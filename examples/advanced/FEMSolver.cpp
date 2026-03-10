/*
 * =====================================================================================
 *
 *       Filename:  FEMSolver.cpp
 *
 *    Description:  FEM Poisson solver for arbitrary 2D meshes:
 *                  triangles (P1), quads (Q1), polygons (centroid-condensed P1),
 *                  embedded flat or on a sphere (Laplace-Beltrami).
 *                  Geometry auto-detected from vertex positions.
 *                  MMS verification: sin(pi x)sin(pi y) for flat,
 *                                    Y1^0 = z/R for sphere.
 *                  MOAB-native data management:
 *                    GLOBAL_ID tag  — 1-based vertex global identifier
 *                    LOCAL_ID  tag  — 0-based DOF index (DofManagerCGP1)
 *                    IS_DIRICHLET / DIRICHLET_VALUE tags on vertices
 *                    DIRICHLET_NODES meshset (ghost-exchange ready)
 *                    Solution / ExactSolution / PointwiseError vertex tags
 *                    ELEM_L2_ERROR element tag for error visualization
 *
 *        Version:  3.0
 *       Compiler:  g++ -std=c++14
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *        Company:  Argonne National Lab
 *
 * =====================================================================================
 */

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cctype>
#include <limits>
#include <cmath>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"

// ─────────────────────────────────────────────────────────────────────────────
namespace util {
    inline void mb_check(moab::ErrorCode ec, const char* msg) {
        if (ec != moab::MB_SUCCESS)
            throw std::runtime_error(std::string("MOAB error: ") + msg
                                     + ", code=" + std::to_string(static_cast<int>(ec)));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
namespace fem {

// =============================================================================
//  Surface Jacobian  (unified for flat-2D and 3D-embedded elements)
// =============================================================================
// For a mapping F: (r,s) -> x in R^3,  J = [F_r | F_s]  is 3x2.
// G = J^T J  (2x2 metric),  sqrt(det G) = area element per unit reference area.
// Surface gradient: grad_S phi_i = J * Ginv * (d phi_i / d xi).
struct SurfJac {
    Eigen::Matrix<double,3,2> J;
    Eigen::Matrix2d            Ginv;
    double                     sqrtDetG;
};

inline SurfJac make_surf_jac(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    SurfJac jac;
    jac.J.col(0) = a;
    jac.J.col(1) = b;
    Eigen::Matrix2d G;
    G(0,0) = a.dot(a);  G(0,1) = a.dot(b);
    G(1,0) = G(0,1);    G(1,1) = b.dot(b);
    // det(G) = |a|²|b|² - (a·b)² = |a×b|².
    // Use the cross-product form: numerically stable even when R >> 1
    // (the algebraic form G₀₀·G₁₁ - G₀₁² suffers catastrophic cancellation
    //  for large-coordinate meshes such as Earth-radius spheres).
    const double detG = a.cross(b).squaredNorm();
    if (detG <= 0.0) throw std::runtime_error("Degenerate element: non-positive metric determinant");
    jac.sqrtDetG = std::sqrt(detG);
    jac.Ginv(0,0) =  G(1,1)/detG;  jac.Ginv(0,1) = -G(0,1)/detG;
    jac.Ginv(1,0) = -G(1,0)/detG;  jac.Ginv(1,1) =  G(0,0)/detG;
    return jac;
}

// Surface gradient of a basis function: 3D tangent vector
inline Eigen::Vector3d surf_grad(const SurfJac& jac, const Eigen::Vector2d& dphi_ref) {
    return jac.J * (jac.Ginv * dphi_ref);
}

// =============================================================================
//  Basis functions and quadrature rules
// =============================================================================

// P1 triangle on reference triangle { (0,0), (1,0), (0,1) }
struct TriP1Basis {
    static constexpr int num_nodes = 3;

    static Eigen::Vector3d eval_phi(const Eigen::Vector2d& xi) {
        return {1.0 - xi[0] - xi[1], xi[0], xi[1]};
    }
    // Returns 3x2: row i = (d phi_i/dr, d phi_i/ds)  — constant
    static Eigen::Matrix<double,3,2> eval_dphi_ref() {
        Eigen::Matrix<double,3,2> d;
        d << -1.0, -1.0,
              1.0,  0.0,
              0.0,  1.0;
        return d;
    }
    // 3-point degree-2 Gauss rule; weights are for the unit-area version,
    // so multiply by 0.5 (reference triangle area) when computing dV.
    struct GaussPoint { Eigen::Vector2d xi; double w; };
    static std::array<GaussPoint,3> gauss_rule() {
        return {{
            {Eigen::Vector2d(1.0/6, 1.0/6), 1.0/3},
            {Eigen::Vector2d(2.0/3, 1.0/6), 1.0/3},
            {Eigen::Vector2d(1.0/6, 2.0/3), 1.0/3}
        }};
    }
};

// Q1 bilinear quad on reference square [-1,1]^2
// Node ordering (CCW): 0=(-1,-1), 1=(+1,-1), 2=(+1,+1), 3=(-1,+1)
struct QuadQ1Basis {
    static constexpr int num_nodes = 4;

    static Eigen::Matrix<double,4,1> eval_phi(double r, double s) {
        Eigen::Matrix<double,4,1> phi;
        phi(0) = 0.25*(1-r)*(1-s);
        phi(1) = 0.25*(1+r)*(1-s);
        phi(2) = 0.25*(1+r)*(1+s);
        phi(3) = 0.25*(1-r)*(1+s);
        return phi;
    }
    // Returns 4x2: row i = (d phi_i/dr, d phi_i/ds) at (r,s)
    static Eigen::Matrix<double,4,2> eval_dphi_ref(double r, double s) {
        Eigen::Matrix<double,4,2> d;
        d(0,0) = -0.25*(1-s);  d(0,1) = -0.25*(1-r);
        d(1,0) =  0.25*(1-s);  d(1,1) = -0.25*(1+r);
        d(2,0) =  0.25*(1+s);  d(2,1) =  0.25*(1+r);
        d(3,0) = -0.25*(1+s);  d(3,1) =  0.25*(1-r);
        return d;
    }
    // 2x2 Gauss rule on [-1,1]^2; weight = 1.0 each (4 points total)
    struct GaussPoint { double r, s, w; };
    static std::array<GaussPoint,4> gauss_rule() {
        static const double g = 1.0 / std::sqrt(3.0);
        return {{
            {-g, -g, 1.0}, { g, -g, 1.0},
            { g,  g, 1.0}, {-g,  g, 1.0}
        }};
    }
};

// =============================================================================
//  File utilities
// =============================================================================
inline bool file_readable(const std::string& p) {
    std::ifstream f(p, std::ios::binary);
    return f.good();
}
inline std::string ext_lower(const std::string& p) {
    auto pos = p.find_last_of('.');
    if (pos == std::string::npos) return {};
    std::string e = p.substr(pos);
    for (char& c : e) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
    return e;
}
inline std::string basename_no_ext(const std::string& p) {
    auto sl = p.find_last_of("/\\");
    std::string n = (sl == std::string::npos) ? p : p.substr(sl+1);
    auto dot = n.find_last_of('.');
    if (dot != std::string::npos) n = n.substr(0, dot);
    return n;
}

// =============================================================================
//  MoabMeshView — wraps moab::Core for any 2D element mesh
// =============================================================================
class MoabMeshView {
public:
    explicit MoabMeshView(const std::string& path) {
        if (!file_readable(path))
            throw std::invalid_argument("Mesh file not readable: " + path);
        const auto ext = ext_lower(path);
        if (ext != ".h5m" && ext != ".vtk" && ext != ".vtu"
                          && ext != ".exo" && ext != ".g")
            throw std::invalid_argument("Unsupported mesh extension: " + ext);

        core_ = std::make_unique<moab::Core>();
        util::mb_check(core_->load_mesh(path.c_str()), "load_mesh");
        setup();
    }

    // Construct from a pre-built Core (used by mesh refinement)
    explicit MoabMeshView(std::unique_ptr<moab::Core> existing_core)
        : core_(std::move(existing_core)) { setup(); }

    const moab::Range& elements() const { return elems_; }
    const moab::Range& vertices()  const { return verts_; }
    const moab::Range& edges()     const { return edges_; }
    moab::Tag           global_id_tag() const { return gid_tag_; }

    // Returns the GLOBAL_ID for a single vertex (1-based, as MOAB convention)
    int global_id(moab::EntityHandle vh) const {
        int gid = -1;
        core_->tag_get_data(gid_tag_, &vh, 1, &gid);
        return gid;
    }

    moab::EntityType elem_type(moab::EntityHandle e) const {
        return core_->type_from_handle(e);
    }

    std::vector<moab::EntityHandle> elem_connectivity(moab::EntityHandle e) const {
        const moab::EntityHandle* ptr = nullptr;
        int n = 0;
        util::mb_check(core_->get_connectivity(e, ptr, n), "elem_connectivity");
        return {ptr, ptr + n};
    }

    std::vector<Eigen::Vector3d> elem_coords(moab::EntityHandle e) const {
        const auto conn = elem_connectivity(e);
        std::vector<double> xyz(3 * conn.size());
        util::mb_check(core_->get_coords(conn.data(), static_cast<int>(conn.size()),
                                          xyz.data()), "elem_coords");
        std::vector<Eigen::Vector3d> out(conn.size());
        for (std::size_t i = 0; i < conn.size(); ++i)
            out[i] = {xyz[3*i], xyz[3*i+1], xyz[3*i+2]};
        return out;
    }

    // Returns edges with exactly one adjacent 2D element (physical boundary)
    moab::Range boundary_edges() const {
        moab::Range bnd;
        for (moab::EntityHandle e : edges_) {
            moab::Range adj;
            core_->get_adjacencies(&e, 1, 2, false, adj);
            if (adj.size() == 1) bnd.insert(e);
        }
        return bnd;
    }

    moab::Core& core() const { return *core_; }

private:
    void setup() {
        util::mb_check(core_->get_entities_by_dimension(0, 0, verts_), "get_verts");
        util::mb_check(core_->get_entities_by_dimension(0, 2, elems_), "get_elems");
        // Create edge adjacencies if missing, then retrieve all edges
        util::mb_check(core_->get_adjacencies(elems_, 1, true, edges_,
                                               moab::Interface::UNION), "get_edges");
        setup_global_ids();
    }

    // Ensure every vertex has a GLOBAL_ID tag.  HDF5 meshes already have one;
    // for VTK and other formats MOAB may not set it, so we assign 1-based
    // sequential IDs ordered by entity handle.
    void setup_global_ids() {
        const int def = 0;
        moab::ErrorCode rc = core_->tag_get_handle(
            "GLOBAL_ID", 1, moab::MB_TYPE_INTEGER, gid_tag_);
        if (rc != moab::MB_SUCCESS || gid_tag_ == 0) {
            // Tag absent: create it
            util::mb_check(
                core_->tag_get_handle("GLOBAL_ID", 1, moab::MB_TYPE_INTEGER, gid_tag_,
                                      moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &def),
                "create GLOBAL_ID tag");
        }
        // Check if vertices already have IDs (HDF5 case)
        bool needs_init = false;
        if (!verts_.empty()) {
            int sample = 0;
            core_->tag_get_data(gid_tag_, &*verts_.begin(), 1, &sample);
            needs_init = (sample == 0);
        }
        if (needs_init) {
            int idx = 1; // 1-based
            for (moab::EntityHandle vh : verts_) {
                util::mb_check(core_->tag_set_data(gid_tag_, &vh, 1, &idx),
                               "set GLOBAL_ID");
                ++idx;
            }
        }
    }

    std::unique_ptr<moab::Core> core_;
    moab::Range verts_, elems_, edges_;
    moab::Tag   gid_tag_ = 0;
};

// =============================================================================
//  Geometry auto-detection
// =============================================================================
enum class GeometryType { FLAT_2D, SPHERE, GENERAL_MANIFOLD };

inline const char* gtype_name(GeometryType g) {
    switch (g) {
        case GeometryType::FLAT_2D:          return "flat 2D";
        case GeometryType::SPHERE:           return "sphere";
        case GeometryType::GENERAL_MANIFOLD: return "general manifold";
    }
    return "unknown";
}

struct GeomInfo { GeometryType type; double sphere_radius; };

static GeomInfo detect_geometry(const MoabMeshView& mesh) {
    const moab::Range& verts = mesh.vertices();
    if (verts.empty()) return {GeometryType::FLAT_2D, 0.0};

    std::vector<double> xyz(3 * verts.size());
    util::mb_check(mesh.core().get_coords(verts, xyz.data()), "detect_geometry coords");

    double max_absz = 0.0, mean_r = 0.0;
    double xmin = xyz[0], xmax = xyz[0];
    double ymin = xyz[1], ymax = xyz[1];
    double zmin = xyz[2], zmax = xyz[2];

    for (std::size_t i = 0; i < verts.size(); ++i) {
        const double x = xyz[3*i], y = xyz[3*i+1], z = xyz[3*i+2];
        max_absz = std::max(max_absz, std::abs(z));
        mean_r  += std::sqrt(x*x + y*y + z*z);
        xmin = std::min(xmin,x); xmax = std::max(xmax,x);
        ymin = std::min(ymin,y); ymax = std::max(ymax,y);
        zmin = std::min(zmin,z); zmax = std::max(zmax,z);
    }
    mean_r /= static_cast<double>(verts.size());

    // Flat if all vertices lie in a plane parallel to xy (z-variation << xy-extent)
    const double xy_extent = std::sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin));
    if ((zmax - zmin) < 1e-6 * std::max(xy_extent, 1e-10))
        return {GeometryType::FLAT_2D, 0.0};

    if (mean_r > 1e-10) {
        double max_rel = 0.0;
        for (std::size_t i = 0; i < verts.size(); ++i) {
            const double x = xyz[3*i], y = xyz[3*i+1], z = xyz[3*i+2];
            max_rel = std::max(max_rel, std::abs(std::sqrt(x*x+y*y+z*z) - mean_r) / mean_r);
        }
        if (max_rel < 1e-3)
            return {GeometryType::SPHERE, mean_r};
    }
    return {GeometryType::GENERAL_MANIFOLD, 0.0};
}

// =============================================================================
//  Unit-sphere rescaling
//  For spherical meshes with R >> 1 (e.g. Earth-radius), normalise all vertex
//  coordinates to the unit sphere before assembly.  This avoids:
//    • catastrophic cancellation in the metric determinant for large R
//    • ill-conditioning of the stiffness matrix
//  The MMS (u = z/R) is the same on the unit sphere (R = 1, u = z).
// =============================================================================
static std::unique_ptr<MoabMeshView>
rescale_sphere(const MoabMeshView& mesh, double R) {
    const moab::Range& verts = mesh.vertices();
    std::vector<double> xyz(3 * verts.size());
    util::mb_check(mesh.core().get_coords(verts, xyz.data()), "rescale:get_coords");
    for (auto& v : xyz) v /= R;

    // Build vertex-handle → index map for connectivity re-wiring
    std::unordered_map<moab::EntityHandle, int> vh_to_idx;
    {
        int idx = 0;
        for (moab::EntityHandle v : verts) vh_to_idx[v] = idx++;
    }

    auto new_core = std::make_unique<moab::Core>();
    moab::Range new_verts;
    util::mb_check(new_core->create_vertices(xyz.data(),
                                              static_cast<int>(verts.size()),
                                              new_verts), "rescale:create_verts");
    std::vector<moab::EntityHandle> vh(new_verts.begin(), new_verts.end());

    for (moab::EntityHandle e : mesh.elements()) {
        const auto conn = mesh.elem_connectivity(e);
        const auto et   = mesh.elem_type(e);
        std::vector<moab::EntityHandle> new_conn(conn.size());
        for (std::size_t k = 0; k < conn.size(); ++k)
            new_conn[k] = vh[vh_to_idx.at(conn[k])];
        moab::EntityHandle new_e;
        util::mb_check(new_core->create_element(et, new_conn.data(),
                                                 static_cast<int>(new_conn.size()),
                                                 new_e), "rescale:create_element");
    }
    return std::make_unique<MoabMeshView>(std::move(new_core));
}

// =============================================================================
//  Mesh refinement — native midpoint subdivision (works for all element types)
//
//  Triangles  → 4 sub-triangles  (red refinement, edge midpoints)
//  Quads      → 4 sub-quads      (edge midpoints + face centroid)
//  Polygons   → n sub-quads      (edge midpoints + face centroid)
//
//  Sphere geometry: new vertices are projected back onto the sphere of radius R.
// =============================================================================
static std::unique_ptr<MoabMeshView>
refine_one_level(const MoabMeshView& mesh, const GeomInfo& ginfo) {
    // --- 1. Extract existing vertex coordinates -------------------------
    const moab::Range& verts = mesh.vertices();
    std::vector<double> xyz_raw(3 * verts.size());
    util::mb_check(mesh.core().get_coords(verts, xyz_raw.data()), "refine:get_coords");

    std::unordered_map<moab::EntityHandle, int> vh_to_idx;
    std::vector<std::array<double, 3>> coords;
    coords.reserve(verts.size());
    {
        int idx = 0;
        for (moab::EntityHandle v : verts) {
            vh_to_idx[v] = idx;
            coords.push_back({xyz_raw[3*idx], xyz_raw[3*idx+1], xyz_raw[3*idx+2]});
            ++idx;
        }
    }

    // --- 2. Edge midpoint cache ----------------------------------------
    // Key: (min_idx, max_idx) → index into coords[]
    std::map<std::pair<int,int>, int> edge_mid_map;

    const bool on_sphere = (ginfo.type == GeometryType::SPHERE);
    const double R = ginfo.sphere_radius;

    auto project_sphere = [&](std::array<double,3>& p) {
        double r = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        if (r > 1e-15) { p[0] *= R/r; p[1] *= R/r; p[2] *= R/r; }
    };

    auto get_or_create_midpoint = [&](int a, int b) -> int {
        if (a > b) std::swap(a, b);
        auto key = std::make_pair(a, b);
        auto it  = edge_mid_map.find(key);
        if (it != edge_mid_map.end()) return it->second;
        std::array<double,3> mp = {
            0.5*(coords[a][0]+coords[b][0]),
            0.5*(coords[a][1]+coords[b][1]),
            0.5*(coords[a][2]+coords[b][2])
        };
        if (on_sphere) project_sphere(mp);
        int new_idx = static_cast<int>(coords.size());
        coords.push_back(mp);
        edge_mid_map[key] = new_idx;
        return new_idx;
    };

    // --- 3. Build refined element list ---------------------------------
    // Reserve ~4× elements (upper bound for tris/quads; polygons may produce more)
    std::vector<std::pair<moab::EntityType, std::vector<int>>> new_elements;
    new_elements.reserve(4 * mesh.elements().size());

    for (moab::EntityHandle e : mesh.elements()) {
        const auto conn = mesh.elem_connectivity(e);
        const int  n    = static_cast<int>(conn.size());
        const auto et   = mesh.elem_type(e);

        std::vector<int> vi(n);
        for (int k = 0; k < n; ++k) vi[k] = vh_to_idx.at(conn[k]);

        // Edge midpoints: mids[k] = midpoint of edge (vi[k], vi[(k+1)%n])
        std::vector<int> mids(n);
        for (int k = 0; k < n; ++k)
            mids[k] = get_or_create_midpoint(vi[k], vi[(k+1)%n]);

        if (et == moab::MBTRI) {
            // Classic red refinement → 4 sub-triangles (all CCW)
            new_elements.push_back({moab::MBTRI, {vi[0],    mids[0], mids[2]}});
            new_elements.push_back({moab::MBTRI, {mids[0],  vi[1],   mids[1]}});
            new_elements.push_back({moab::MBTRI, {mids[2],  mids[1], vi[2]  }});
            new_elements.push_back({moab::MBTRI, {mids[0],  mids[1], mids[2]}});
        } else if (et == moab::MBQUAD) {
            // Quad → 4 sub-quads via face centroid + edge midpoints
            std::array<double,3> cen = {0,0,0};
            for (int k = 0; k < 4; ++k) {
                cen[0] += coords[vi[k]][0];
                cen[1] += coords[vi[k]][1];
                cen[2] += coords[vi[k]][2];
            }
            cen[0] /= 4; cen[1] /= 4; cen[2] /= 4;
            if (on_sphere) project_sphere(cen);
            int ci = static_cast<int>(coords.size());
            coords.push_back(cen);

            new_elements.push_back({moab::MBQUAD, {vi[0],   mids[0], ci,      mids[3]}});
            new_elements.push_back({moab::MBQUAD, {mids[0], vi[1],   mids[1], ci     }});
            new_elements.push_back({moab::MBQUAD, {ci,      mids[1], vi[2],   mids[2]}});
            new_elements.push_back({moab::MBQUAD, {mids[3], ci,      mids[2], vi[3]  }});
        } else {
            // MBPOLYGON (n-gon) → n sub-quads via face centroid + edge midpoints
            // Each sub-quad k: {vi[k], mids[k], centroid, mids[(k-1+n)%n]}  (CCW)
            std::array<double,3> cen = {0,0,0};
            for (int k = 0; k < n; ++k) {
                cen[0] += coords[vi[k]][0];
                cen[1] += coords[vi[k]][1];
                cen[2] += coords[vi[k]][2];
            }
            cen[0] /= n; cen[1] /= n; cen[2] /= n;
            if (on_sphere) project_sphere(cen);
            int ci = static_cast<int>(coords.size());
            coords.push_back(cen);

            for (int k = 0; k < n; ++k)
                new_elements.push_back({moab::MBQUAD,
                    {vi[k], mids[k], ci, mids[(k-1+n)%n]}});
        }
    }

    // --- 4. Build new MOAB Core ----------------------------------------
    auto new_core = std::make_unique<moab::Core>();

    std::vector<double> flat_coords;
    flat_coords.reserve(3 * coords.size());
    for (const auto& c : coords) {
        flat_coords.push_back(c[0]);
        flat_coords.push_back(c[1]);
        flat_coords.push_back(c[2]);
    }

    moab::Range new_verts;
    util::mb_check(new_core->create_vertices(flat_coords.data(),
                                              static_cast<int>(coords.size()),
                                              new_verts), "refine:create_vertices");

    // Convert Range to vector for O(1) indexed access
    std::vector<moab::EntityHandle> vh(new_verts.begin(), new_verts.end());

    for (std::size_t ei = 0; ei < new_elements.size(); ++ei) {
        const moab::EntityType             et   = new_elements[ei].first;
        const std::vector<int>&            conn = new_elements[ei].second;
        std::vector<moab::EntityHandle> ch(conn.size());
        for (std::size_t k = 0; k < conn.size(); ++k) ch[k] = vh[conn[k]];
        moab::EntityHandle new_e;
        util::mb_check(new_core->create_element(et, ch.data(),
                                                 static_cast<int>(ch.size()), new_e),
                       "refine:create_element");
    }

    return std::make_unique<MoabMeshView>(std::move(new_core));
}

// =============================================================================
//  Manufactured Solutions (polymorphic)
// =============================================================================
struct MMSBase {
    virtual double          u     (double x, double y, double z) const = 0;
    virtual double          rhs_f (double x, double y, double z) const = 0;
    // Exact gradient in 3D; for flat MMS z-component is zero
    virtual Eigen::Vector3d grad  (double x, double y, double z) const = 0;
    virtual ~MMSBase() = default;
};

// Flat 2D: -kappa * Delta u = f,  u = sin(pi x) sin(pi y)
struct MMSinSin : MMSBase {
    double kappa = 1.0;
    double u(double x, double y, double /*z*/) const override {
        return std::sin(M_PI*x) * std::sin(M_PI*y);
    }
    Eigen::Vector3d grad(double x, double y, double /*z*/) const override {
        return { M_PI*std::cos(M_PI*x)*std::sin(M_PI*y),
                 M_PI*std::sin(M_PI*x)*std::cos(M_PI*y),
                 0.0 };
    }
    double rhs_f(double x, double y, double z) const override {
        return kappa * 2.0 * M_PI * M_PI * u(x, y, z);
    }
};

// Sphere: -Delta_S u = f,  u = z/R  (spherical harmonic Y1^0)
//   Delta_S u = -(2/R^2) u   =>   f = (2/R^2) * u = 2z/R^3
struct MMSSphere : MMSBase {
    double R;
    explicit MMSSphere(double r) : R(r) {}
    double u(double /*x*/, double /*y*/, double z) const override {
        return z / R;
    }
    // Surface gradient of z/R on sphere:
    //   grad_S(z/R) = (1/R)(e_z - (z/R) n_hat)  where n_hat = (x,y,z)/R
    Eigen::Vector3d grad(double x, double y, double z) const override {
        const double R2 = R*R;
        return { -x*z/(R2*R), -y*z/(R2*R), (R2 - z*z)/(R2*R) };
    }
    double rhs_f(double /*x*/, double /*y*/, double z) const override {
        return 2.0 * z / (R*R*R);
    }
};

// =============================================================================
//  DOF manager (vertex-based CG, works for any element type)
// =============================================================================
// =============================================================================
//  DOF manager  —  CG P1: one DOF per vertex
//  Uses a MOAB dense integer tag "LOCAL_ID" (0-based) on vertices so that
//  DOF indices live natively in the mesh database.  This is the foundation
//  for future MPI-parallel assembly: the same tag can be ghost-exchanged via
//  MOAB's parallel communicator to give remote processors their ghost DOF
//  indices without any additional bookkeeping.
// =============================================================================
class DofManagerCGP1 {
public:
    explicit DofManagerCGP1(const MoabMeshView& mesh) : mesh_(mesh) {
        moab::Core& core = mesh_.core();
        const int def = -1;
        // Create (or retrieve) a dense integer tag "LOCAL_ID" on all vertices.
        // Using MB_TAG_DENSE ensures every vertex always has a value; default -1
        // makes uninitialized entries detectable.
        util::mb_check(
            core.tag_get_handle("LOCAL_ID", 1, moab::MB_TYPE_INTEGER, local_id_tag_,
                                moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &def),
            "create LOCAL_ID tag");

        // Assign sequential 0-based DOF indices in Range order (deterministic,
        // matches entity-handle ordering which is stable for a given Core).
        int idx = 0;
        for (moab::EntityHandle vh : mesh_.vertices()) {
            util::mb_check(core.tag_set_data(local_id_tag_, &vh, 1, &idx),
                           "set LOCAL_ID");
            ++idx;
        }
        ndofs_ = idx;
    }

    int      ndofs()         const { return ndofs_; }
    moab::Tag local_id_tag() const { return local_id_tag_; }

    // Single-vertex DOF lookup — reads tag from MOAB Core.
    int vertex_dof(moab::EntityHandle vh) const {
        int dof = -1;
        util::mb_check(mesh_.core().tag_get_data(local_id_tag_, &vh, 1, &dof),
                       "vertex_dof");
        return dof;
    }

    // Element DOFs — batch tag read on all connectivity handles in one call,
    // which is significantly faster than N individual tag lookups for large
    // meshes and avoids hash-map overhead entirely.
    std::vector<int> elem_dofs(moab::EntityHandle e) const {
        const auto conn = mesh_.elem_connectivity(e);
        std::vector<int> d(conn.size(), -1);
        util::mb_check(
            mesh_.core().tag_get_data(local_id_tag_, conn.data(),
                                      static_cast<int>(conn.size()), d.data()),
            "elem_dofs");
        return d;
    }

private:
    const MoabMeshView& mesh_;
    moab::Tag           local_id_tag_ = 0;
    int                 ndofs_{};
};

// =============================================================================
//  Strong Dirichlet enforcement with consistent RHS correction
// =============================================================================
static void apply_dirichlet_strong(Eigen::SparseMatrix<double>& A,
                                   Eigen::VectorXd& b,
                                   std::vector<char>& is_bc,
                                   std::vector<double>& gvals) {
    const int N = static_cast<int>(is_bc.size());
    int bc_count = 0;
    for (int i = 0; i < N; ++i) bc_count += is_bc[i] ? 1 : 0;
    if (bc_count == 0) {
        is_bc[0] = 1;
        gvals[0] = b[0];
        std::cerr << "Warning: no Dirichlet DOFs; anchoring DOF 0.\n";
    }

    A.makeCompressed();

    // Pass 1: zero BC columns, subtract A(i,j)*g_j from b[i], set diagonal to 1
    for (int j = 0; j < A.outerSize(); ++j) {
        if (!is_bc[j]) continue;
        bool diag_seen = false;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A,j); it; ++it) {
            if (it.row() != j) {
                b[it.row()] -= it.value() * gvals[j];
                it.valueRef() = 0.0;
            } else {
                it.valueRef() = 1.0;
                diag_seen = true;
            }
        }
        if (!diag_seen) {
            // Diagonal missing: rebuild with it inserted
            std::vector<Eigen::Triplet<double>> trips;
            trips.reserve(A.nonZeros() + 1);
            for (int c = 0; c < A.outerSize(); ++c)
                for (Eigen::SparseMatrix<double>::InnerIterator it(A,c); it; ++it)
                    trips.emplace_back(it.row(), it.col(), it.value());
            trips.emplace_back(j, j, 1.0);
            Eigen::SparseMatrix<double> Anew(N,N);
            Anew.setFromTriplets(trips.begin(), trips.end());
            A.swap(Anew);
        }
    }

    // Pass 2: zero off-diagonal entries in BC rows
    for (int j = 0; j < A.outerSize(); ++j) {
        if (is_bc[j]) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A,j); it; ++it)
            if (is_bc[it.row()] && it.row() != j)
                it.valueRef() = 0.0;
    }

    for (int i = 0; i < N; ++i)
        if (is_bc[i]) b[i] = gvals[i];
}

// =============================================================================
//  Poisson assembler  —  kappa * int grad_S phi_i . grad_S phi_j dS = int f phi_i dS
// =============================================================================
class PoissonAssembler {
public:
    PoissonAssembler(const MoabMeshView& m, const DofManagerCGP1& d, double k,
                     const GeomInfo& g = GeomInfo{GeometryType::FLAT_2D, 1.0})
        : mesh_(m), dofs_(d), kappa_(k), ginfo_(g) {}

    void assemble(Eigen::SparseMatrix<double>& A,
                  Eigen::VectorXd& b,
                  const MMSBase& mms) const {
        const int N = dofs_.ndofs();
        std::vector<Eigen::Triplet<double>> trips;
        b = Eigen::VectorXd::Zero(N);

        for (moab::EntityHandle e : mesh_.elements()) {
            switch (mesh_.elem_type(e)) {
                case moab::MBTRI:     assemble_tri    (e, mms, trips, b); break;
                case moab::MBQUAD:    assemble_quad   (e, mms, trips, b); break;
                case moab::MBPOLYGON: assemble_polygon(e, mms, trips, b); break;
                default: break; // skip 3D entities if any
            }
        }
        A.resize(N,N);
        A.setFromTriplets(trips.begin(), trips.end());
    }

private:
    const MoabMeshView&   mesh_;
    const DofManagerCGP1& dofs_;
    double                kappa_;
    GeomInfo              ginfo_;

    // --- P1 triangle ---
    void assemble_tri(moab::EntityHandle e,
                      const MMSBase& mms,
                      std::vector<Eigen::Triplet<double>>& trips,
                      Eigen::VectorXd& b) const {
        const auto X     = mesh_.elem_coords(e);   // 3 nodes
        const auto gdofs = dofs_.elem_dofs(e);
        const SurfJac jac = make_surf_jac(X[1]-X[0], X[2]-X[0]);
        const auto dphi   = TriP1Basis::eval_dphi_ref();

        Eigen::Matrix3d Ke = Eigen::Matrix3d::Zero();
        Eigen::Vector3d fe = Eigen::Vector3d::Zero();

        for (const auto& gp : TriP1Basis::gauss_rule()) {
            const double dV  = gp.w * 0.5 * jac.sqrtDetG;
            const Eigen::Vector3d xq = X[0] + jac.J * gp.xi;
            const double fval = mms.rhs_f(xq[0], xq[1], xq[2]);
            const Eigen::Vector3d phi = TriP1Basis::eval_phi(gp.xi);
            for (int i = 0; i < 3; ++i) {
                const Eigen::Vector3d gi = surf_grad(jac, dphi.row(i).transpose());
                fe[i] += phi[i] * fval * dV;
                for (int j = 0; j < 3; ++j)
                    Ke(i,j) += kappa_ * gi.dot(surf_grad(jac, dphi.row(j).transpose())) * dV;
            }
        }
        scatter3(Ke, fe, gdofs, trips, b);
    }

    // --- Q1 quad ---
    void assemble_quad(moab::EntityHandle e,
                       const MMSBase& mms,
                       std::vector<Eigen::Triplet<double>>& trips,
                       Eigen::VectorXd& b) const {
        const auto X     = mesh_.elem_coords(e);   // 4 nodes
        const auto gdofs = dofs_.elem_dofs(e);

        Eigen::Matrix<double,4,4> Ke = Eigen::Matrix<double,4,4>::Zero();
        Eigen::Matrix<double,4,1> fe = Eigen::Matrix<double,4,1>::Zero();

        for (const auto& gp : QuadQ1Basis::gauss_rule()) {
            const auto dphi = QuadQ1Basis::eval_dphi_ref(gp.r, gp.s);
            const auto phi  = QuadQ1Basis::eval_phi(gp.r, gp.s);
            // Isoparametric Jacobian and physical coordinate at this Gauss point
            Eigen::Vector3d Jcol0 = Eigen::Vector3d::Zero();
            Eigen::Vector3d Jcol1 = Eigen::Vector3d::Zero();
            Eigen::Vector3d xq    = Eigen::Vector3d::Zero();
            for (int i = 0; i < 4; ++i) {
                Jcol0 += X[i] * dphi(i,0);
                Jcol1 += X[i] * dphi(i,1);
                xq    += X[i] * phi(i);
            }
            const SurfJac jac = make_surf_jac(Jcol0, Jcol1);
            const double dV   = gp.w * jac.sqrtDetG;
            const double fval = mms.rhs_f(xq[0], xq[1], xq[2]);
            for (int i = 0; i < 4; ++i) {
                const Eigen::Vector3d gi = surf_grad(jac, dphi.row(i).transpose());
                fe[i] += phi(i) * fval * dV;
                for (int j = 0; j < 4; ++j)
                    Ke(i,j) += kappa_ * gi.dot(surf_grad(jac, dphi.row(j).transpose())) * dV;
            }
        }
        scatter4(Ke, fe, gdofs, trips, b);
    }

    // --- Polygon: centroid-based triangulation with static condensation ---
    // Each n-gon is split into n sub-triangles {centroid, v_k, v_{k+1}}.
    // The centroid DOF is locally condensed out before scattering to the global
    // system, so no global centroid DOF is introduced.  This avoids the degenerate
    // sub-triangles that can arise from fan-triangulation (v0 as origin) when three
    // nearly-collinear vertices appear across the polygon.
    void assemble_polygon(moab::EntityHandle e,
                          const MMSBase& mms,
                          std::vector<Eigen::Triplet<double>>& trips,
                          Eigen::VectorXd& b) const {
        const auto X     = mesh_.elem_coords(e);
        const auto gdofs = dofs_.elem_dofs(e);
        const int  n     = static_cast<int>(X.size());
        if (n < 3) return;

        // Centroid: mean of vertices, projected onto sphere if needed
        Eigen::Vector3d cpos = Eigen::Vector3d::Zero();
        for (const auto& xi : X) cpos += xi;
        cpos /= static_cast<double>(n);
        if (ginfo_.type == GeometryType::SPHERE) {
            const double r = cpos.norm();
            if (r > 1e-12) cpos *= ginfo_.sphere_radius / r;
        }

        // Build (n+1)×(n+1) local system.  Local DOF 0 = centroid; 1..n = polygon vertices.
        const int m = n + 1;
        Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(m, m);
        Eigen::VectorXd fe = Eigen::VectorXd::Zero(m);

        const auto dphi = TriP1Basis::eval_dphi_ref();

        for (int k = 0; k < n; ++k) {
            const int kp1 = (k + 1) % n;
            const std::array<Eigen::Vector3d,3> Xs = {cpos, X[k], X[kp1]};
            const std::array<int,3>             ls = {0, k+1, kp1+1};
            const SurfJac jac = make_surf_jac(Xs[1]-Xs[0], Xs[2]-Xs[0]);
            if (jac.sqrtDetG < 1e-14) continue;

            for (const auto& gp : TriP1Basis::gauss_rule()) {
                const double dV  = gp.w * 0.5 * jac.sqrtDetG;
                const Eigen::Vector3d xq   = Xs[0] + jac.J * gp.xi;
                const double fval = mms.rhs_f(xq[0], xq[1], xq[2]);
                const Eigen::Vector3d phi  = TriP1Basis::eval_phi(gp.xi);
                for (int i = 0; i < 3; ++i) {
                    const Eigen::Vector3d gi = surf_grad(jac, dphi.row(i).transpose());
                    fe[ls[i]] += phi[i] * fval * dV;
                    for (int j = 0; j < 3; ++j)
                        Ke(ls[i], ls[j]) += kappa_ * gi.dot(surf_grad(jac, dphi.row(j).transpose())) * dV;
                }
            }
        }

        // Static condensation: eliminate centroid DOF (local index 0).
        // K_red[i,j] = Ke[i,j] - Ke[i,0]*Ke[0,j]/Ke[0,0]   (i,j = 1..n)
        // f_red[i]   = fe[i]   - Ke[i,0]*fe[0]  /Ke[0,0]
        const double Kcc = Ke(0, 0);
        if (std::abs(Kcc) < 1e-30) return;

        for (int i = 1; i < m; ++i) {
            fe[i] -= Ke(i,0) * fe[0] / Kcc;
            for (int j = 1; j < m; ++j)
                Ke(i,j) -= Ke(i,0) * Ke(0,j) / Kcc;
        }

        // Scatter condensed n×n system to global
        for (int i = 1; i < m; ++i) {
            b[gdofs[i-1]] += fe[i];
            for (int j = 1; j < m; ++j)
                trips.emplace_back(gdofs[i-1], gdofs[j-1], Ke(i,j));
        }
    }

    void scatter3(const Eigen::Matrix3d& Ke, const Eigen::Vector3d& fe,
                  const std::vector<int>& d,
                  std::vector<Eigen::Triplet<double>>& trips,
                  Eigen::VectorXd& b) const {
        for (int i = 0; i < 3; ++i) {
            b[d[i]] += fe[i];
            for (int j = 0; j < 3; ++j) trips.emplace_back(d[i],d[j],Ke(i,j));
        }
    }
    void scatter4(const Eigen::Matrix<double,4,4>& Ke,
                  const Eigen::Matrix<double,4,1>& fe,
                  const std::vector<int>& d,
                  std::vector<Eigen::Triplet<double>>& trips,
                  Eigen::VectorXd& b) const {
        for (int i = 0; i < 4; ++i) {
            b[d[i]] += fe[i];
            for (int j = 0; j < 4; ++j) trips.emplace_back(d[i],d[j],Ke(i,j));
        }
    }
};

// =============================================================================
//  Dirichlet BC setup — geometry-aware
//  Persists BC information in two MOAB tags on vertices:
//    IS_DIRICHLET   (int,    dense, default 0)  — 1 for Dirichlet nodes
//    DIRICHLET_VALUE(double, dense, default 0.0) — prescribed value
//  Also creates (or replaces) a meshset "DIRICHLET_NODES" containing exactly
//  the Dirichlet vertices, ready for ghost exchange in a future parallel run.
//  The out-parameters is_bc / gvals are LOCAL_ID-indexed arrays for immediate
//  use in apply_dirichlet_strong(); they are derived from the tags.
// =============================================================================
static void build_dirichlet(const MoabMeshView& mesh,
                             const DofManagerCGP1& dofs,
                             const MMSBase& mms,
                             GeometryType gtype,
                             std::vector<char>& is_bc,
                             std::vector<double>& gvals) {
    moab::Core& core = mesh.core();
    const int N = dofs.ndofs();

    // ── Create (or get) IS_DIRICHLET and DIRICHLET_VALUE tags ──────────────
    moab::Tag is_dir_tag, dir_val_tag;
    const int    def_int = 0;
    const double def_dbl = 0.0;
    util::mb_check(
        core.tag_get_handle("IS_DIRICHLET", 1, moab::MB_TYPE_INTEGER, is_dir_tag,
                            moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &def_int),
        "IS_DIRICHLET tag");
    util::mb_check(
        core.tag_get_handle("DIRICHLET_VALUE", 1, moab::MB_TYPE_DOUBLE, dir_val_tag,
                            moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &def_dbl),
        "DIRICHLET_VALUE tag");

    // Reset tags to defaults on all vertices (ensures a clean state at each level)
    {
        std::vector<moab::EntityHandle> vhs(mesh.vertices().begin(),
                                             mesh.vertices().end());
        std::vector<int>    zeros(vhs.size(), 0);
        std::vector<double> zerods(vhs.size(), 0.0);
        core.tag_set_data(is_dir_tag,  vhs.data(), static_cast<int>(vhs.size()), zeros.data());
        core.tag_set_data(dir_val_tag, vhs.data(), static_cast<int>(vhs.size()), zerods.data());
    }

    // ── Helper: mark one vertex as Dirichlet via tags ───────────────────────
    auto mark_bc_vertex = [&](moab::EntityHandle vh, double val) {
        const int one = 1;
        core.tag_set_data(is_dir_tag,  &vh, 1, &one);
        core.tag_set_data(dir_val_tag, &vh, 1, &val);
    };

    // ── Determine which vertices are Dirichlet ───────────────────────────────
    if (gtype == GeometryType::SPHERE) {
        const moab::Range bnd = mesh.boundary_edges();
        if (bnd.empty()) {
            // Closed sphere: pin one vertex to remove Laplace-Beltrami null space
            moab::EntityHandle vh = *mesh.vertices().begin();
            double xyz[3] = {};
            util::mb_check(core.get_coords(&vh, 1, xyz), "sphere anchor coords");
            mark_bc_vertex(vh, mms.u(xyz[0], xyz[1], xyz[2]));
            std::cout << "  [BC] Sphere: pinned DOF " << dofs.vertex_dof(vh)
                      << " to remove Laplace-Beltrami null space.\n";
        } else {
            // Open (regional) sphere: Dirichlet on all physical boundary vertices
            std::cout << "  [BC] Sphere with boundary: applying Dirichlet on "
                      << bnd.size() << " boundary edges.\n";
            for (moab::EntityHandle e : bnd) {
                const moab::EntityHandle* vhs = nullptr; int nconn = 0;
                util::mb_check(core.get_connectivity(e, vhs, nconn), "edge conn");
                if (nconn != 2 || vhs == nullptr) continue;
                double xyz[6] = {};
                util::mb_check(core.get_coords(vhs, 2, xyz), "edge vcoords");
                for (int k = 0; k < 2; ++k) {
                    int flag = 0;
                    core.tag_get_data(is_dir_tag, &vhs[k], 1, &flag);
                    if (!flag)
                        mark_bc_vertex(vhs[k], mms.u(xyz[3*k], xyz[3*k+1], xyz[3*k+2]));
                }
            }
        }
    } else {
        // Flat 2D or general manifold: Dirichlet on all physical boundary vertices
        const moab::Range bnd = mesh.boundary_edges();
        int marked = 0;
        for (moab::EntityHandle e : bnd) {
            const moab::EntityHandle* vhs = nullptr; int nconn = 0;
            util::mb_check(core.get_connectivity(e, vhs, nconn), "edge conn");
            if (nconn != 2 || vhs == nullptr) continue;
            double xyz[6] = {};
            util::mb_check(core.get_coords(vhs, 2, xyz), "edge vcoords");
            for (int k = 0; k < 2; ++k) {
                int flag = 0;
                core.tag_get_data(is_dir_tag, &vhs[k], 1, &flag);
                if (!flag) {
                    mark_bc_vertex(vhs[k], mms.u(xyz[3*k], xyz[3*k+1], xyz[3*k+2]));
                    ++marked;
                }
            }
        }
        if (marked == 0) {
            moab::EntityHandle vh = *mesh.vertices().begin();
            double xyz[3] = {};
            util::mb_check(core.get_coords(&vh, 1, xyz), "anchor coords");
            mark_bc_vertex(vh, mms.u(xyz[0], xyz[1], xyz[2]));
            std::cerr << "Warning: no boundary detected; anchored DOF "
                      << dofs.vertex_dof(vh) << ".\n";
        }
    }

    // ── Build DIRICHLET_NODES meshset (for ghost exchange in parallel) ───────
    // Delete old set if it exists, then create a fresh one.
    {
        moab::Tag set_tag;
        if (core.tag_get_handle("DIRICHLET_NODES_SET_HANDLE", 1,
                                moab::MB_TYPE_HANDLE, set_tag) == moab::MB_SUCCESS) {
            moab::EntityHandle old_set = 0;
            core.tag_get_data(set_tag, NULL, 0, &old_set);
            if (old_set) core.delete_entities(&old_set, 1);
        }
    }
    moab::EntityHandle bc_meshset;
    util::mb_check(core.create_meshset(moab::MESHSET_SET, bc_meshset),
                   "create DIRICHLET_NODES meshset");
    // Collect Dirichlet vertices by scanning the tag
    for (moab::EntityHandle vh : mesh.vertices()) {
        int flag = 0;
        core.tag_get_data(is_dir_tag, &vh, 1, &flag);
        if (flag) core.add_entities(bc_meshset, &vh, 1);
    }
    // Store the meshset handle in a tag so it can be retrieved later
    moab::Tag set_handle_tag;
    const moab::EntityHandle def_handle = 0;
    util::mb_check(
        core.tag_get_handle("DIRICHLET_NODES_SET_HANDLE", 1, moab::MB_TYPE_HANDLE,
                            set_handle_tag,
                            moab::MB_TAG_MESH | moab::MB_TAG_CREAT, &def_handle),
        "DIRICHLET_NODES_SET_HANDLE tag");
    util::mb_check(core.tag_set_data(set_handle_tag, NULL, 0, &bc_meshset),
                   "store bc_meshset handle");

    // ── Derive LOCAL_ID-indexed arrays for the linear-algebra layer ──────────
    is_bc.assign(static_cast<std::size_t>(N), 0);
    gvals.assign(static_cast<std::size_t>(N), 0.0);
    for (moab::EntityHandle vh : mesh.vertices()) {
        int flag = 0;  double val = 0.0;
        core.tag_get_data(is_dir_tag,  &vh, 1, &flag);
        core.tag_get_data(dir_val_tag, &vh, 1, &val);
        if (flag) {
            const int d = dofs.vertex_dof(vh);
            is_bc [static_cast<std::size_t>(d)] = 1;
            gvals [static_cast<std::size_t>(d)] = val;
        }
    }
}

// =============================================================================
//  Error computation (L2 and H1-seminorm)
//  Also stores the per-element L2-error contribution in a MOAB tag
//  "ELEM_L2_ERROR" (dense double) on the 2D elements, which can be loaded
//  into ParaView to visualize where the solution error is concentrated.
// =============================================================================
struct ErrorMetrics { double L2{}, H1semi{}, h{}; };

static ErrorMetrics compute_errors(const MoabMeshView& mesh,
                                    const DofManagerCGP1& dofs,
                                    const Eigen::VectorXd& x,
                                    const MMSBase& mms,
                                    const GeomInfo& ginfo = GeomInfo{GeometryType::FLAT_2D,1.0}) {
    moab::Core& core = mesh.core();

    // Create (or get) element-wise L2 error tag for visualization
    moab::Tag elem_err_tag;
    const double def_err = 0.0;
    util::mb_check(
        core.tag_get_handle("ELEM_L2_ERROR", 1, moab::MB_TYPE_DOUBLE, elem_err_tag,
                            moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &def_err),
        "ELEM_L2_ERROR tag");

    double L2sq = 0.0, H1sq = 0.0, sum_h = 0.0;
    int count = 0;

    const auto qr_tri  = TriP1Basis::gauss_rule();
    const auto qr_quad = QuadQ1Basis::gauss_rule();
    const auto dphi_tri = TriP1Basis::eval_dphi_ref();

    for (moab::EntityHandle e : mesh.elements()) {
        const moab::EntityType et = mesh.elem_type(e);
        double elem_L2sq = 0.0;

        if (et == moab::MBTRI) {
            const auto X     = mesh.elem_coords(e);
            const auto gdofs = dofs.elem_dofs(e);
            Eigen::Vector3d uh;
            uh << x[gdofs[0]], x[gdofs[1]], x[gdofs[2]];

            const SurfJac jac = make_surf_jac(X[1]-X[0], X[2]-X[0]);
            Eigen::Vector3d grad_uh = Eigen::Vector3d::Zero();
            for (int i = 0; i < 3; ++i)
                grad_uh += uh[i] * surf_grad(jac, dphi_tri.row(i).transpose());

            for (const auto& gp : qr_tri) {
                const double dV  = gp.w * 0.5 * jac.sqrtDetG;
                const Eigen::Vector3d xq = X[0] + jac.J * gp.xi;
                const double uex = mms.u(xq[0], xq[1], xq[2]);
                const double uhq = uh.dot(TriP1Basis::eval_phi(gp.xi));
                elem_L2sq += (uex - uhq)*(uex - uhq) * dV;
            }
            L2sq += elem_L2sq;
            const Eigen::Vector3d xmid = (X[0]+X[1]+X[2]) / 3.0;
            H1sq += (mms.grad(xmid[0],xmid[1],xmid[2]) - grad_uh).squaredNorm()
                     * (0.5 * jac.sqrtDetG);

            sum_h += std::max({(X[1]-X[0]).norm(), (X[2]-X[1]).norm(), (X[0]-X[2]).norm()});
            ++count;

        } else if (et == moab::MBQUAD) {
            const auto X     = mesh.elem_coords(e);
            const auto gdofs = dofs.elem_dofs(e);
            Eigen::Matrix<double,4,1> uh;
            for (int i = 0; i < 4; ++i) uh(i) = x[gdofs[i]];

            for (const auto& gp : qr_quad) {
                const auto dphi = QuadQ1Basis::eval_dphi_ref(gp.r, gp.s);
                const auto phi  = QuadQ1Basis::eval_phi(gp.r, gp.s);
                Eigen::Vector3d Jcol0 = Eigen::Vector3d::Zero();
                Eigen::Vector3d Jcol1 = Eigen::Vector3d::Zero();
                Eigen::Vector3d xq    = Eigen::Vector3d::Zero();
                for (int i = 0; i < 4; ++i) {
                    Jcol0 += X[i] * dphi(i,0);
                    Jcol1 += X[i] * dphi(i,1);
                    xq    += X[i] * phi(i);
                }
                const SurfJac jac = make_surf_jac(Jcol0, Jcol1);
                const double dV   = gp.w * jac.sqrtDetG;
                const double uex  = mms.u(xq[0], xq[1], xq[2]);
                const double uhq  = uh.dot(phi);
                elem_L2sq += (uex - uhq)*(uex - uhq) * dV;
                Eigen::Vector3d grad_uh = Eigen::Vector3d::Zero();
                for (int i = 0; i < 4; ++i)
                    grad_uh += uh(i) * surf_grad(jac, dphi.row(i).transpose());
                H1sq += (mms.grad(xq[0],xq[1],xq[2]) - grad_uh).squaredNorm() * dV;
            }
            L2sq += elem_L2sq;
            sum_h += std::max({(X[1]-X[0]).norm(),(X[2]-X[1]).norm(),
                               (X[3]-X[2]).norm(),(X[0]-X[3]).norm()});
            ++count;

        } else if (et == moab::MBPOLYGON) {
            const auto X     = mesh.elem_coords(e);
            const auto gdofs = dofs.elem_dofs(e);
            const int  n     = static_cast<int>(X.size());
            if (n < 3) continue;

            Eigen::Vector3d cpos = Eigen::Vector3d::Zero();
            for (const auto& xi : X) cpos += xi;
            cpos /= static_cast<double>(n);
            if (ginfo.type == GeometryType::SPHERE) {
                const double r = cpos.norm();
                if (r > 1e-12) cpos *= ginfo.sphere_radius / r;
            }
            double uh_c = 0.0;
            for (int k = 0; k < n; ++k) uh_c += x[gdofs[k]];
            uh_c /= static_cast<double>(n);

            double h_poly    = 0.0;
            double poly_L2sq = 0.0; // accumulate polygon's own L2sq for tag
            for (int k = 0; k < n; ++k) {
                const int kp1 = (k + 1) % n;
                const std::array<Eigen::Vector3d,3> Xs = {cpos, X[k], X[kp1]};
                Eigen::Vector3d uh; uh << uh_c, x[gdofs[k]], x[gdofs[kp1]];

                const SurfJac jac = make_surf_jac(Xs[1]-Xs[0], Xs[2]-Xs[0]);
                if (jac.sqrtDetG < 1e-14) continue;
                Eigen::Vector3d grad_uh = Eigen::Vector3d::Zero();
                for (int i = 0; i < 3; ++i)
                    grad_uh += uh[i] * surf_grad(jac, dphi_tri.row(i).transpose());

                double subtri_L2sq = 0.0;
                for (const auto& gp : qr_tri) {
                    const double dV  = gp.w * 0.5 * jac.sqrtDetG;
                    const Eigen::Vector3d xq = Xs[0] + jac.J * gp.xi;
                    const double uex = mms.u(xq[0], xq[1], xq[2]);
                    const double uhq = uh.dot(TriP1Basis::eval_phi(gp.xi));
                    subtri_L2sq += (uex - uhq)*(uex - uhq) * dV;
                }
                poly_L2sq += subtri_L2sq;
                L2sq      += subtri_L2sq;
                const Eigen::Vector3d xmid = (Xs[0]+Xs[1]+Xs[2]) / 3.0;
                H1sq += (mms.grad(xmid[0],xmid[1],xmid[2]) - grad_uh).squaredNorm()
                         * (0.5 * jac.sqrtDetG);
                h_poly = std::max({h_poly, (Xs[1]-Xs[0]).norm(),
                                   (Xs[2]-Xs[1]).norm(), (Xs[0]-Xs[2]).norm()});
            }
            sum_h += h_poly;
            ++count;
            // Set per-element L2 error tag for polygon elements
            {
                const double eL2 = std::sqrt(poly_L2sq);
                core.tag_set_data(elem_err_tag, &e, 1, &eL2);
            }
            elem_L2sq = 0.0; // reset (unused in polygon path beyond this point)
        }

        // Store per-element L2 error contribution (sqrt) in MOAB tag
        if (et == moab::MBTRI || et == moab::MBQUAD) {
            const double eL2 = std::sqrt(elem_L2sq);
            core.tag_set_data(elem_err_tag, &e, 1, &eL2);
        }
    }
    ErrorMetrics em;
    em.L2     = std::sqrt(L2sq);
    em.H1semi = std::sqrt(H1sq);
    em.h      = (count > 0) ? sum_h / static_cast<double>(count) : 0.0;
    return em;
}

// =============================================================================
//  Convergence rate (log-log slope)
// =============================================================================
static double compute_rate(const std::vector<double>& hs,
                            const std::vector<double>& errs) {
    const int n = static_cast<int>(hs.size());
    if (n < 2) return 0.0;
    double sx=0, sy=0, sxx=0, sxy=0;
    for (int i = 0; i < n; ++i) {
        const double lx = std::log(hs[i]), ly = std::log(errs[i]);
        sx+=lx; sy+=ly; sxx+=lx*lx; sxy+=lx*ly;
    }
    const double d = n*sxx - sx*sx;
    return (std::abs(d) < 1e-14) ? 0.0 : (n*sxy - sx*sy)/d;
}

// =============================================================================
//  Solution tag storage and VTK export
// =============================================================================

// Store the FEM solution, exact MMS solution, and pointwise error as MOAB
// dense tags on every mesh vertex immediately after the linear solve.
// Having this data in MOAB tags means it is:
//   (a) accessible from any downstream code without passing Eigen vectors, and
//   (b) ready for ghost exchange once MOAB parallel comm is in place.
static void store_solution_tags(const MoabMeshView& mesh,
                                 const DofManagerCGP1& dofs,
                                 const Eigen::VectorXd& x,
                                 const MMSBase& mms) {
    moab::Core& core = mesh.core();
    std::vector<moab::EntityHandle> vhs(mesh.vertices().begin(),
                                         mesh.vertices().end());
    const int nv = static_cast<int>(vhs.size());

    std::vector<double> xyz(3 * nv);
    util::mb_check(core.get_coords(vhs.data(), nv, xyz.data()), "store_tags:coords");

    std::vector<double> sol(nv), exact(nv), err(nv);
    for (int i = 0; i < nv; ++i) {
        const double uh  = x[dofs.vertex_dof(vhs[i])];
        const double uex = mms.u(xyz[3*i], xyz[3*i+1], xyz[3*i+2]);
        sol  [i] = uh;
        exact[i] = uex;
        err  [i] = std::abs(uh - uex);
    }

    const double def = 0.0;
    auto set_tag = [&](const char* name, const std::vector<double>& vals) {
        moab::Tag tag;
        util::mb_check(
            core.tag_get_handle(name, 1, moab::MB_TYPE_DOUBLE, tag,
                                moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &def),
            name);
        util::mb_check(core.tag_set_data(tag, vhs.data(), nv, vals.data()), name);
    };
    set_tag("Solution",       sol);
    set_tag("ExactSolution",  exact);
    set_tag("PointwiseError", err);
}

// Write the mesh with all currently-set tags to a VTK file.
// Only the 2D surface elements are exported (edges created for adjacency
// queries are excluded so ParaView sees a proper 3D surface mesh).
static void write_solution(const MoabMeshView& mesh,
                            const std::string& out_base) {
    moab::Core& core = mesh.core();
    const std::string out = out_base + "_solution.vtk";
    moab::EntityHandle eset;
    util::mb_check(core.create_meshset(moab::MESHSET_SET, eset), "create_meshset");
    util::mb_check(core.add_entities(eset, mesh.elements()), "add_entities");
    util::mb_check(core.write_file(out.c_str(), nullptr, nullptr, &eset, 1), "write_file");
    util::mb_check(core.delete_entities(&eset, 1), "delete meshset");
}

} // namespace fem

// =============================================================================
//  main
// =============================================================================
// Shared solve helper: assemble → BCs → direct/iterative solve → error metrics
static fem::ErrorMetrics run_level(const fem::MoabMeshView& mesh,
                                   const fem::GeomInfo&     ginfo,
                                   const fem::MMSBase&      mms,
                                   const std::string&       vtk_base = "") {
    fem::DofManagerCGP1 cgdofs(mesh);
    std::cout << "  DOFs=" << cgdofs.ndofs()
              << "  Elems=" << mesh.elements().size() << "\n";

    fem::PoissonAssembler asmbl(mesh, cgdofs, 1.0, ginfo);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    asmbl.assemble(A, b, mms);

    std::vector<char>   is_bc;
    std::vector<double> gvals;
    fem::build_dirichlet(mesh, cgdofs, mms, ginfo.type, is_bc, gvals);
    fem::apply_dirichlet_strong(A, b, is_bc, gvals);

    Eigen::VectorXd x;
    bool solved = false;

    // Prefer CG with incomplete-Cholesky preconditioner (memory-efficient, scalable).
    // Relative tolerance 1e-14; fall back to direct LLT if CG fails to converge.
    {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower|Eigen::Upper,
                                 Eigen::IncompleteCholesky<double>> cg;
        cg.setTolerance(1e-14);
        cg.setMaxIterations(std::max(10000, 10 * cgdofs.ndofs()));
        cg.compute(A);
        x = cg.solve(b);
        solved = (cg.info() == Eigen::Success);
        if (solved)
            std::cout << "  Solver: CG  iters=" << cg.iterations()
                      << "  res=" << cg.error() << "\n";
        else
            std::cout << "  CG did not converge (res=" << cg.error()
                      << "); falling back to direct LLT.\n";
    }
    if (!solved) {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
        llt.compute(A);
        if (llt.info() == Eigen::Success) {
            x = llt.solve(b);
            solved = (llt.info() == Eigen::Success);
            if (solved) std::cout << "  Solver: direct LLT\n";
        }
    }
    if (!solved) throw std::runtime_error("Linear solver failed");

    // Store solution fields as MOAB tags on mesh vertices immediately after solve.
    // compute_errors also stores ELEM_L2_ERROR on elements.
    fem::store_solution_tags(mesh, cgdofs, x, mms);
    const auto em = fem::compute_errors(mesh, cgdofs, x, mms, ginfo);
    std::cout << "  h~" << em.h
              << "  |L2|=" << em.L2
              << "  |H1-semi|=" << em.H1semi << "\n";

    if (!vtk_base.empty()) {
        fem::write_solution(mesh, vtk_base);
        std::cout << "  Wrote: " << vtk_base << "_solution.vtk\n";
    }
    return em;
}

int main(int argc, char** argv) {
    try {
        // --- Argument parsing -----------------------------------------
        int num_levels = 0;
        std::vector<std::string> mesh_paths;

        for (int ai = 1; ai < argc; ++ai) {
            std::string arg = argv[ai];
            if ((arg == "-l" || arg == "--levels") && ai + 1 < argc) {
                num_levels = std::stoi(argv[++ai]);
            } else {
                mesh_paths.push_back(arg);
            }
        }

        if (mesh_paths.empty()) {
            std::cerr << "Usage: " << argv[0]
                      << " [-l NUM_LEVELS] mesh1 [mesh2 ...]\n"
                      << "  -l N : Auto-refine mesh N times for convergence study\n"
                      << "  (No -l) : supply meshes coarse→fine manually\n";
            return 1;
        }

        std::vector<double> hs, L2s, H1s;

        if (num_levels > 0) {
            // === Refinement-hierarchy mode ================================
            const std::string base_path = mesh_paths[0];
            std::cout << "\n=== Refinement study: " << base_path
                      << "  (" << num_levels << " levels) ===\n";

            fem::MoabMeshView base_mesh(base_path);
            auto ginfo = fem::detect_geometry(base_mesh);
            std::cout << "  Geometry : " << fem::gtype_name(ginfo.type);
            if (ginfo.type == fem::GeometryType::SPHERE)
                std::cout << " (R=" << ginfo.sphere_radius << ")";
            std::cout << "\n";

            // Rescale spherical meshes to unit sphere for numerical stability
            std::unique_ptr<fem::MoabMeshView> rescaled_base;
            const fem::MoabMeshView* working_base = &base_mesh;
            if (ginfo.type == fem::GeometryType::SPHERE
                    && std::abs(ginfo.sphere_radius - 1.0) > 1e-3) {
                std::cout << "  Rescaling to unit sphere (R="
                          << ginfo.sphere_radius << " → 1)\n";
                rescaled_base = fem::rescale_sphere(base_mesh, ginfo.sphere_radius);
                ginfo.sphere_radius = 1.0;
                working_base = rescaled_base.get();
            }

            // Build the refinement hierarchy: level 0 = (rescaled) original mesh
            std::vector<const fem::MoabMeshView*> levels;
            std::vector<std::unique_ptr<fem::MoabMeshView>> refined;
            levels.push_back(working_base);
            for (int lv = 0; lv < num_levels; ++lv) {
                std::cout << "  Refining level " << lv << " → " << lv+1 << " ...\n";
                refined.push_back(fem::refine_one_level(*levels.back(), ginfo));
                levels.push_back(refined.back().get());
            }

            // MMS fixed for all levels (geometry already known)
            std::unique_ptr<fem::MMSBase> mms;
            if (ginfo.type == fem::GeometryType::SPHERE) {
                mms = std::make_unique<fem::MMSSphere>(ginfo.sphere_radius);
                std::cout << "  MMS : u = z/R  (Y1^0 spherical harmonic)\n";
            } else {
                auto fm = std::make_unique<fem::MMSinSin>(); fm->kappa = 1.0;
                mms = std::move(fm);
                std::cout << "  MMS : u = sin(pi x) sin(pi y)\n";
            }

            // Solve on each level, export VTK only for the finest
            for (int lv = 0; lv <= num_levels; ++lv) {
                std::cout << "\n  --- Level " << lv << " ---\n";
                const bool is_finest = (lv == num_levels);
                const std::string vtk = is_finest
                    ? fem::basename_no_ext(base_path) + "_level" + std::to_string(lv)
                    : "";
                const auto em = run_level(*levels[lv], ginfo, *mms, vtk);
                hs .push_back(em.h);
                L2s.push_back(em.L2);
                H1s.push_back(em.H1semi);
            }

        } else {
            // === Manual mesh-list mode (original behaviour) ==============
            for (int ai = 0; ai < static_cast<int>(mesh_paths.size()); ++ai) {
                const std::string& path = mesh_paths[ai];
                std::cout << "\n=== Mesh: " << path << " ===\n";

                fem::MoabMeshView mesh(path);
                auto ginfo = fem::detect_geometry(mesh);
                std::cout << "  Geometry : " << fem::gtype_name(ginfo.type);
                if (ginfo.type == fem::GeometryType::SPHERE)
                    std::cout << "  (R=" << ginfo.sphere_radius << ")";
                std::cout << "\n";

                // Rescale to unit sphere for numerical stability
                std::unique_ptr<fem::MoabMeshView> rescaled;
                const fem::MoabMeshView* working_mesh = &mesh;
                if (ginfo.type == fem::GeometryType::SPHERE
                        && std::abs(ginfo.sphere_radius - 1.0) > 1e-3) {
                    std::cout << "  Rescaling to unit sphere (R="
                              << ginfo.sphere_radius << " → 1)\n";
                    rescaled = fem::rescale_sphere(mesh, ginfo.sphere_radius);
                    ginfo.sphere_radius = 1.0;
                    working_mesh = rescaled.get();
                }

                std::unique_ptr<fem::MMSBase> mms;
                if (ginfo.type == fem::GeometryType::SPHERE) {
                    mms = std::make_unique<fem::MMSSphere>(ginfo.sphere_radius);
                    std::cout << "  MMS : u = z/R  (Y1^0)\n";
                } else {
                    auto fm = std::make_unique<fem::MMSinSin>(); fm->kappa = 1.0;
                    mms = std::move(fm);
                    std::cout << "  MMS : u = sin(pi x) sin(pi y)\n";
                }

                const bool is_last = (ai == static_cast<int>(mesh_paths.size()) - 1);
                const std::string vtk = is_last ? fem::basename_no_ext(path) : "";
                const auto em = run_level(*working_mesh, ginfo, *mms, vtk);
                hs .push_back(em.h);
                L2s.push_back(em.L2);
                H1s.push_back(em.H1semi);
            }
        }

        if (hs.size() >= 2) {
            std::cout << "\nConvergence (log-log):"
                      << "  L2 rate="      << fem::compute_rate(hs, L2s)
                      << "  H1-semi rate=" << fem::compute_rate(hs, H1s) << "\n";
        } else {
            std::cout << "\nProvide multiple meshes or use -l N for convergence rates.\n";
        }

    } catch (const std::exception& ex) {
        std::cerr << "Fatal: " << ex.what() << "\n";
        return 2;
    }
    return 0;
}

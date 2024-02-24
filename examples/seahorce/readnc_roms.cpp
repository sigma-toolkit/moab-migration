#include <iostream>
#include <vector>
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include <netcdfcpp.h>

#define get_2d_tag_data(tag_name, tag_data) \
    { \
        moab::Tag fieldTag; \
        result = twodMesh.tag_get_handle(tag_name, 1, moab::MB_TYPE_DOUBLE, fieldTag, moab::MB_TAG_DENSE);MB_CHK_ERR(result); \
        tag_data.resize(ents2d.size());\
        result = twodMesh.tag_get_data(fieldTag, ents2d, tag_data.data());MB_CHK_ERR(result); \
    }

#define get_3d_tag_data(tag_name, tag_data) \
    { \
        moab::Tag fieldTag; \
        result = threedMesh.tag_get_handle(tag_name, 1, moab::MB_TYPE_DOUBLE, fieldTag, moab::MB_TAG_DENSE);MB_CHK_ERR(result); \
        tag_data.resize(ents3d.size());\
        result = threedMesh.tag_get_data(fieldTag, ents3d, tag_data.data());MB_CHK_ERR(result); \
    }


int main() {
    moab::ErrorCode result;
    // Start MOAB instances
    moab::Core twodMesh;
    moab::Core threedMesh;

    // Create meshsets
    moab::EntityHandle twod_set, threed_set;

    // Create meshsets
    result = twodMesh.create_meshset(moab::MESHSET_SET, twod_set);MB_CHK_ERR(result);
    result = threedMesh.create_meshset(moab::MESHSET_SET, threed_set);MB_CHK_ERR(result);

    // Load the files
    result = twodMesh.load_file("roms_2d_projected.h5m", &twod_set);MB_CHK_ERR(result);
    result = threedMesh.load_file("roms_3d_projected.h5m", &threed_set);MB_CHK_ERR(result);

    // Query the root set for all 2D and 3D entities
    moab::Range verts2d, ents2d, ents3d;
    result = twodMesh.get_entities_by_dimension(twod_set, 0, verts2d);MB_CHK_ERR(result);
    result = twodMesh.get_entities_by_dimension(twod_set, 2, ents2d);MB_CHK_ERR(result);
    result = threedMesh.get_entities_by_dimension(threed_set, 3, ents3d);MB_CHK_ERR(result);

    // Get the bathymetry tag
    std::vector<double> bathymetry, ssh;
    get_2d_tag_data("Bathymetry", bathymetry);
    get_2d_tag_data("SSH", ssh);

    // Get the velocity tag
    moab::Tag velocityy_tag;
    result = threedMesh.tag_get_handle("ROMS_VelMeridional", 1, moab::MB_TYPE_DOUBLE, velocityy_tag, moab::MB_TAG_DENSE);MB_CHK_ERR(result);

    // Get the velocity data
    std::vector<double> velocityy(ents3d.size());
    result = threedMesh.tag_get_data(velocityy_tag, ents3d, &velocityy[0]);MB_CHK_ERR(result);

    // TODO: Reshape and average the velocity data

    // Create a new NetCDF file
    NcFile ncfile("roms_his.nc", NcFile::Replace);

    // Set the title
    const char* titlestr="Wind-Driven Upwelling/Downwelling for the North-Atlantic test case";
    constexpr int xi_rho_init = 413;
    constexpr int eta_rho_init = 147;
    constexpr int xi_u_init = 412;
    constexpr int eta_u_init = 147;
    constexpr int xi_v_init = 413;
    constexpr int eta_v_init = 146;
    constexpr int s_rho_init = 100;
    // s_rho_init = int(ents3d.size()/ents2d.size())

    // Create dimensions
    NcDim* xi_rho = ncfile.add_dim("xi_rho", xi_rho_init);
    NcDim* eta_rho = ncfile.add_dim("eta_rho", eta_rho_init);
    NcDim* s_rho = ncfile.add_dim("s_rho", s_rho_init);
    NcDim* xi_u = ncfile.add_dim("xi_u", xi_u_init);
    NcDim* eta_u = ncfile.add_dim("eta_u", eta_u_init);
    NcDim* xi_v = ncfile.add_dim("xi_v", xi_v_init);
    NcDim* eta_v = ncfile.add_dim("eta_v", eta_v_init);
    NcDim* ocean_time = ncfile.add_dim("ocean_time");

    // Set attributes
    ncfile.add_att("title", titlestr);
    ncfile.add_att("type", "ROMS/TOMS restart file");

    // Create variables
    NcVar* theta_s = ncfile.add_var("theta_s", ncDouble);
    theta_s->add_att("long_name", "S-coordinate surface control parameter");
    theta_s->put(std::vector<double>{5.0}.data());

    NcVar* theta_b = ncfile.add_var("theta_b", ncDouble);
    theta_b->add_att("long_name", "S-coordinate bottom control parameter");
    theta_b->put(std::vector<double>{0.5}.data());

    NcVar* hc = ncfile.add_var("hc", ncDouble);
    hc->add_att("long_name", "S-coordinate parameter, critical depth");
    hc->add_att("units", "meter");
    hc->put(std::vector<double>{100.0}.data());

    NcVar* bath_var = ncfile.add_var("h", ncDouble, eta_rho, xi_rho);
    bath_var->add_att("long_name", "bathymetry at RHO-points");
    bath_var->add_att("units", "meter");
    bath_var->add_att("grid", "grid");
    bath_var->add_att("location", "face");
    bath_var->add_att("coordinates", "lon_rho lat_rho");
    bath_var->add_att("field", "bath, scalar");
    std::cout << "Bathymetry values: " << bathymetry.size() << ": " << bathymetry[0] << ", " << bathymetry[1] << std::endl;
    bath_var->put(bathymetry.data(), eta_rho_init, xi_rho_init);

    NcVar* zeta_var = ncfile.add_var("zeta", ncDouble, ocean_time, eta_rho, xi_rho);
    zeta_var->add_att("long_name", "free-surface");
    zeta_var->add_att("units", "meter");
    zeta_var->add_att("time", "ocean_time");
    zeta_var->add_att("grid", "grid");
    zeta_var->add_att("location", "face");
    zeta_var->add_att("coordinates", "x_rho y_rho ocean_time");
    zeta_var->add_att("field", "free-surface, scalar, series");
    std::cout << "Zeta values: " << ssh.size() << ": " << ssh[0] << ", " << ssh[1] << std::endl;
    zeta_var->put(ssh.data(), 1, eta_rho_init, xi_rho_init);

    ncfile.close();

    return 0;
}

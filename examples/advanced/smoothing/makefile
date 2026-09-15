# Include the MOAB configuration information so that
# all required flags and libs are populated correctly
include ../../makefile.config

default: all

# LaplacianSmoother calls MPI_Init unconditionally and reads/writes the mesh in
# parallel, so it needs both MPI and HDF5.
ALLEXAMPLES =
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
ALLEXAMPLES += LaplacianSmoother
endif

# LloydRelaxation builds against MPI alone but reads an .h5m and writes
# lloydfinal.h5m, so only its run rule needs HDF5.  This is the gating it had in
# examples/advanced/makefile, kept as-is by the move.
ifeq ("$(MOAB_MPI_ENABLED)","yes")
ALLEXAMPLES += LloydRelaxation
endif

# The optimization examples take their gradients from Eigen's unsupported
# AutoDiff module, so Eigen3 is not optional for them.  They generate their own
# meshes and so need no HDF5, but they do call MPI_Init unconditionally.
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_EIGEN3_ENABLED)","yes-yes")
ALLEXAMPLES += MeshOptimization SphericalCVTOptimization
endif

all: $(ALLEXAMPLES)

LaplacianSmoother: LaplacianSmoother.o ${MOAB_LIBFILE}
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-LaplacianSmoother: LaplacianSmoother
	${RUNPARALLEL} ./LaplacianSmoother -n 5

LloydRelaxation: LloydRelaxation.o ${MOAB_LIBFILE}
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-LloydRelaxation: LloydRelaxation
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./LloydRelaxation
	${RUNPARALLEL} ./LloydRelaxation
endif

MeshOptimization: MeshOptimization.o ${MOAB_LIBFILE}
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

# The self-test is the regression check: closed-form metric values, the Verdict
# cross-check and a finite-difference check of the AutoDiff gradient.
run-MeshOptimization: MeshOptimization
	${RUNSERIAL} ./MeshOptimization --selftest
	${RUNPARALLEL} ./MeshOptimization -d 2 -N 8
	${RUNPARALLEL} ./MeshOptimization -d 3 -N 4

SphericalCVTOptimization: SphericalCVTOptimization.o ${MOAB_LIBFILE}
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-SphericalCVTOptimization: SphericalCVTOptimization
	${RUNPARALLEL} ./SphericalCVTOptimization -N 200 -l 20 -n 30

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

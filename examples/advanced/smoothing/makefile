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

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

# OptimizeMeshMesquite is deliberately not built.  It drives Mesquite through
# the iMesh bindings, and neither is present: there is no itaps/ directory in
# this tree and no build system here detects Mesquite.  The source is kept for
# reference until it is either ported or removed.

all: $(ALLEXAMPLES)

LaplacianSmoother: LaplacianSmoother.o ${MOAB_LIBFILE}
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ $< $(MOAB_LIBS_LINK)

run-LaplacianSmoother: LaplacianSmoother
	${RUNPARALLEL} ./LaplacianSmoother -n 5

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

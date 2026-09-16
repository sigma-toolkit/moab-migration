# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

# Both examples call the iMOAB Fortran module, which lives in libMOAB itself.
# DirectAccessF90 builds its mesh in memory and needs no MPI; MigrateMesh is
# parallel by nature.
F90EXAMPLES = DirectAccessF90
F90PAREXAMPLES = MigrateMesh

ALLEXAMPLES =
ifeq ("$(MOAB_FORTRAN_ENABLED)","yes")
ALLEXAMPLES += ${F90EXAMPLES}
endif

ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_FORTRAN_ENABLED)","yes-yes")
parallel: ${F90PAREXAMPLES}
ALLEXAMPLES += ${F90PAREXAMPLES}
else
parallel:
endif

all: $(ALLEXAMPLES)

DirectAccessF90: DirectAccessF90.o ${MOAB_LIBFILE}
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${MOAB_LIBS_LINK} ${MOAB_CXX_RUNTIME_LIBS} ${MOAB_FC_MPI_LIBS}

run-DirectAccessF90: DirectAccessF90
	${RUNSERIAL} ./DirectAccessF90

MigrateMesh: MigrateMesh.o ${MOAB_LIBFILE}
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${MOAB_LIBS_LINK} ${MOAB_CXX_RUNTIME_LIBS} ${MOAB_FC_MPI_LIBS}

run-MigrateMesh: MigrateMesh
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./MigrateMesh
endif

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

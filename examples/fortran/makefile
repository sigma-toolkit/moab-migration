# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

# DirectAccessNoHolesF90 and PushParMeshIntoMoabF90 go through iMesh
# (iMesh_f.h, iMeshP_f.h) and link ${IMESH_LIBS}.  DirectAccessF90 and
# MigrateMesh only call the iMOAB Fortran module, which lives in libMOAB
# itself, so they must not be gated on iMesh - doing so was why this directory
# built nothing at all.  DirectAccessF90 additionally needs no MPI.
F90EXAMPLES = DirectAccessNoHolesF90
F90PAREXAMPLES = PushParMeshIntoMoabF90
IMOAB_F90EXAMPLES = DirectAccessF90
IMOAB_F90PAREXAMPLES = MigrateMesh

ALLEXAMPLES =
ifeq ("$(MOAB_FORTRAN_ENABLED)","yes")
ALLEXAMPLES += ${IMOAB_F90EXAMPLES}
endif

ifeq ("$(MOAB_FORTRAN_ENABLED)-$(MOAB_IMESH_ENABLED)","yes-yes")
ALLEXAMPLES += ${F90EXAMPLES}
endif

ifeq ("$(MOAB_MPI_ENABLED)","yes")
ifeq ("$(MOAB_FORTRAN_ENABLED)","yes")
parallel: ${IMOAB_F90PAREXAMPLES}
ALLEXAMPLES += ${IMOAB_F90PAREXAMPLES}
endif
ifeq ("$(MOAB_FORTRAN_ENABLED)-$(MOAB_IMESH_ENABLED)","yes-yes")
parallel: ${F90PAREXAMPLES}
ALLEXAMPLES += ${F90PAREXAMPLES}
endif
else
parallel:
endif

all: $(ALLEXAMPLES)

DirectAccessNoHolesF90: DirectAccessNoHolesF90.o ${MOAB_LIBFILE}
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${IMESH_LIBS}

run-DirectAccessNoHolesF90: DirectAccessNoHolesF90
	${RUNSERIAL} ./DirectAccessNoHolesF90 -n 1000

PushParMeshIntoMoabF90: PushParMeshIntoMoabF90.o
	@echo "[F90LD]  $@"
	${VERBOSE}${MOAB_FC} -o $@ $< ${IMESH_LIBS}

run-PushParMeshIntoMoabF90: PushParMeshIntoMoabF90
	${RUNSERIAL} ./PushParMeshIntoMoabF90
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./PushParMeshIntoMoabF90
endif

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


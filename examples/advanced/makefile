# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES =  ExtrudePoly GeneralGnomonic AddEdgeFluxVector
PAREXAMPLES = LloydRelaxation CrystalRouterExample GenLargeMesh
EXOIIEXAMPLES = VisTags
ALLEXAMPLES = ${EXAMPLES}

ifeq ("$(MOAB_NETCDF_ENABLED)","yes")
ALLEXAMPLES += $(EXOIIEXAMPLES)
endif

ifeq ("$(MOAB_MPI_ENABLED)","yes")
parallel: ${PAREXAMPLES}
ALLEXAMPLES += ${PAREXAMPLES}
ifeq ("$(ENABLE_mbcoupler)","yes")
PAREXAMPLES += CoupleMGen
endif
else
parallel:
endif

all: $(ALLEXAMPLES)

ExtrudePoly: ExtrudePoly.o ${MOAB_LIBFILE}
LloydRelaxation: LloydRelaxation.o ${MOAB_LIBFILE}
CrystalRouterExample: CrystalRouterExample.o ${MOAB_LIBFILE}
VisTags: VisTags.o ${MOAB_LIBFILE}
GenLargeMesh: GenLargeMesh.o ${MOAB_LIBFILE}
GeneralGnomonic: GeneralGnomonic.o ${MOAB_LIBFILE}
AddEdgeFluxVector: AddEdgeFluxVector.o ${MOAB_LIBFILE}

ExtrudePoly LloydRelaxation CrystalRouterExample VisTags GenLargeMesh GeneralGnomonic AddEdgeFluxVector:
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

CoupleMGen: CoupleMGen.o ${MOAB_LIBFILE}
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< -lmbcoupler ${MOAB_LIBS_LINK}

run-ExtrudePoly: ExtrudePoly
	${RUNSERIAL} ./ExtrudePoly

run-LloydRelaxation: LloydRelaxation
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./LloydRelaxation
	${RUNPARALLEL} ./LloydRelaxation
endif

run-CrystalRouterExample: CrystalRouterExample
	${RUNSERIAL} ./CrystalRouterExample
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./CrystalRouterExample -n 2 -t 10 -r 0
endif

run-VisTags: VisTags
	${RUNSERIAL} ./VisTags

run-GenLargeMesh: GenLargeMesh
	${RUNSERIAL} ./GenLargeMesh
	${RUNSERIAL} ./GenLargeMesh -A 2 -B 2 -C 5
	${RUNSERIAL} ./GenLargeMesh -A 2 -B 2 -C 5 -f
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./GenLargeMesh -M $(NPROCS) -N 1 -K 1 -A 1 -B 2 -C 2 -r
	${RUNPARALLEL} ./GenLargeMesh -M $(NPROCS) -N 1 -K 1 -A 1 -B 2 -C 2 -d 10.0 -b 2 -r
	${RUNPARALLEL} ./GenLargeMesh -M $(NPROCS) -N 1 -K 1 -A 1 -B 2 -C 2 -d 10.0 -b 2 -r -f
	${RUNPARALLEL} ./GenLargeMesh -M $(NPROCS)  -b 10  -A 1 -C 1 -B 1 -r -G
endif

run-CoupleMGen: CoupleMGen
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./CoupleMGen -M $(NPROCS)
endif


run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}


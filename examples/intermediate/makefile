# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES = DirectAccessWithHoles DirectAccessNoHoles PointInElementSearch ExtractClose QuadTriConv FixPolys
PAREXAMPLES = ReduceExchangeTags ExchangeHalos

ifeq ("$(MOAB_HDF5_ENABLED)","yes")
EXAMPLES += LoadPartial ComputeTriDual
endif
ALLEXAMPLES = ${EXAMPLES}

ifeq ("$(MOAB_MPI_ENABLED)","yes")
parallel: ${PAREXAMPLES}
ALLEXAMPLES += ${PAREXAMPLES}
else
parallel:
endif

all: $(ALLEXAMPLES)

LoadPartial: LoadPartial.o ${MOAB_LIBFILE}
QuadTriConv: QuadTriConv.o ${MOAB_LIBFILE}
DirectAccessWithHoles: DirectAccessWithHoles.o ${MOAB_LIBFILE}
DirectAccessNoHoles: DirectAccessNoHoles.o ${MOAB_LIBFILE}
ReduceExchangeTags: ReduceExchangeTags.o ${MOAB_LIBFILE}
PointInElementSearch: PointInElementSearch.o ${MOAB_LIBFILE}
ExtractClose: ExtractClose.o ${MOAB_LIBFILE}
FixPolys: FixPolys.o ${MOAB_LIBFILE}
ComputeTriDual: ComputeTriDual.o ${MOAB_LIBFILE}
ExchangeHalos: ExchangeHalos.o ${MOAB_LIBFILE}

$(EXAMPLES) $(PAREXAMPLES):
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-LoadPartial: LoadPartial
	${RUNSERIAL} ./LoadPartial

run-DirectAccessWithHoles: DirectAccessWithHoles
	${RUNSERIAL} ./DirectAccessWithHoles -n 1000 -H 5

run-DirectAccessNoHoles: DirectAccessNoHoles
	${RUNSERIAL} ./DirectAccessNoHoles -n 1000

run-ReduceExchangeTags: ReduceExchangeTags
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL} ./ReduceExchangeTags
	${RUNPARALLEL} ./ReduceExchangeTags
endif

run-PointInElementSearch: PointInElementSearch
	${RUNSERIAL} ./PointInElementSearch
ifeq ("$(MOAB_MPI_ENABLED)","yes")
	${RUNPARALLEL} ./PointInElementSearch
endif

run-ExtractClose: ExtractClose 
	${RUNSERIAL} ./ExtractClose -h

run-FixPolys: FixPolys
	${RUNSERIAL} ./FixPolys -h

run-QuadTriConv: QuadTriConv
	${RUNSERIAL} ./QuadTriConv

run-ComputeTriDual: ComputeTriDual
	${RUNSERIAL} ./ComputeTriDual

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

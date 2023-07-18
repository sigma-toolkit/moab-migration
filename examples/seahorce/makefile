# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES = 
PAREXAMPLES = 
#MOAB_CXXFLAGS += -Xpreprocessor -fopenmp #-g3 -O0
MOAB_CXXFLAGS += -g3

ifeq ("$(MOAB_MPI_ENABLED)","yes")
ifeq ("$(MOAB_HDF5_ENABLED)","yes")
PAREXAMPLES += RemapMPASROMS FilterElements
parallel: ${PAREXAMPLES}
endif
else
parallel:
endif

ALLEXAMPLES = ${EXAMPLES} ${PAREXAMPLES}
all: $(ALLEXAMPLES)

RemapMPASROMS: RemapMPASROMS.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

FilterElements: FilterElements.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

#$(EXAMPLES) $(PAREXAMPLES):
#	@echo "[CXXLD]  $@"
#	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

clean: clobber
	rm -rf ${ALLEXAMPLES}

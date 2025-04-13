# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	MOAB_CXXFLAGS += -fopenmp
	MOAB_LIBS_LINK += -fopenmp
endif
#MOAB_CXXFLAGS += -std=c++17 -Wno-missing-declarations -DKOKKOS_DEPENDENCE -Iinclude  -I/opt/compadre/include -I/opt/compadre/include/kokkos -Iexternal/alglib-cpp/src/
#MOAB_CXXFLAGS += -g3 -O0 -Iinclude
MOAB_CXXFLAGS += -g -O3 -Iinclude

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. 
PREDCXXFLAGS = -O0
EXTLIBS = lib/libtet.a # external/alglib-cpp/src/libalglib.a
#COMPADRELIBS=-L/opt/compadre/lib -lcompadre -lkokkoskernels -lkokkoscore

EXAMPLES = RemapMPASROMS FilterElements
all: RemapMPASROMS

RemapMPASROMS.o: include/*.hpp
RemapMPASROMS: RemapMPASROMS.o ${MOAB_LIBDIR}/libMOAB.la ${EXTLIBS}
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} ${EXTLIBS}

readnc_roms.o: include/*.hpp
readnc_roms: readnc_roms.o ${MOAB_LIBDIR}/libMOAB.la ${EXTLIBS}
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} ${EXTLIBS}

lib/libTSFIT.a: 
	${VERBOSE}make -C src/tsfit

predicates.o: src/predicates.cxx
	${VERBOSE}${MOAB_CXX} $(PREDCXXFLAGS) -Iinclude -c src/predicates.cxx

FilterElements: FilterElements.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

Test1DInterpolant: Test1DInterpolant.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestDelaunay: TestDelaunay.o ${MOAB_LIBDIR}/libMOAB.la  lib/libtet.a
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}  lib/libtet.a

define compile-command
${VERBOSE}${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<
endef

%.o: src/%.cpp 
	@echo "  [CXX]  $<"
	$(compile-command)

clean: clobber
	rm -rf ${EXAMPLES}

# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	MOAB_CXXFLAGS += -fopenmp
	MOAB_LIBS_LINK += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
	MOAB_CXXFLAGS += -Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/17.0.3/include
	MOAB_LIBS_LINK += -L/opt/homebrew/Cellar/libomp/17.0.3/lib -lomp
endif
MOAB_CXXFLAGS += -std=c++17 -Wno-missing-declarations -DKOKKOS_DEPENDENCE -Iinclude -I/opt/compadre/include -I/opt/compadre/include/kokkos -Iexternal/alglib-cpp/src/
#MOAB_CXXFLAGS += -g3 -O0

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. 
PREDCXXFLAGS = -O0

EXAMPLES = RemapMPASROMS FilterElements
all: RemapMPASROMS

RemapMPASROMS.o: include/*.hpp
RemapMPASROMS: RemapMPASROMS.o ${MOAB_LIBDIR}/libMOAB.la lib/libtet.a lib/libTSFIT.a
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} lib/libtet.a  lib/libTSFIT.a -L/opt/compadre/lib -lcompadre -lkokkoskernels -lkokkoscore external/alglib-cpp/src/libalglib.a

lib/libtet.a: src/tetgen.cxx predicates.o
	${VERBOSE}${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} -c src/tetgen.cxx
	ar r lib/libtet.a tetgen.o predicates.o

lib/libTSFIT.a: src/tetgen.cxx predicates.o
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

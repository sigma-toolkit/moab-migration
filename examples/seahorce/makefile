# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	MOAB_CXXFLAGS += -fopenmp
	MOAB_LIBS_LINK += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
	MOAB_CXXFLAGS += -Xpreprocessor -fopenmp
	MOAB_LIBS_LINK += -L/opt/homebrew/Cellar/libomp/16.0.6/lib -lomp
endif
MOAB_CXXFLAGS += -Iinclude

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. 
PREDCXXFLAGS = -O0

EXAMPLES = RemapMPASROMS FilterElements
all: RemapMPASROMS

RemapMPASROMS.o: include/*.hpp
RemapMPASROMS: RemapMPASROMS.o ${MOAB_LIBDIR}/libMOAB.la libtet.a
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}  -fopenmp libtet.a

libtet.a: src/tetgen.cxx predicates.o
	${VERBOSE}${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} -DTETLIBRARY -c src/tetgen.cxx
	ar r libtet.a tetgen.o predicates.o

predicates.o: src/predicates.cxx
	${VERBOSE}${MOAB_CXX} $(PREDCXXFLAGS) -Iinclude -c src/predicates.cxx

FilterElements: FilterElements.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

Test1DInterpolant: Test1DInterpolant.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestDelaunay: TestDelaunay.o ${MOAB_LIBDIR}/libMOAB.la  libtet.a
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}  libtet.a -fopenmp

define compile-command
${VERBOSE}${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<
endef

%.o: src/%.cpp 
	@echo "  [CXX]  $<"
	$(compile-command)

clean: clobber
	rm -rf ${EXAMPLES}

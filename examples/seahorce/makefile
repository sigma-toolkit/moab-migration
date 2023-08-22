# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

MOAB_CXXFLAGS += -Xpreprocessor -fopenmp #-g3 -O0
MOAB_CXXFLAGS += -g3 -Iinclude

EXAMPLES = RemapMPASROMS FilterElements
all: RemapMPASROMS

RemapMPASROMS: RemapMPASROMS.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} -L/opt/homebrew/Cellar/libomp/16.0.6/lib -lomp

FilterElements: FilterElements.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

pchip_test: pchip_test.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< 

define compile-command
${VERBOSE}${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<
endef

%.o: src/%.cpp
	@echo "  [CXX]  $<"
	$(compile-command)

clean: clobber
	rm -rf ${EXAMPLES}

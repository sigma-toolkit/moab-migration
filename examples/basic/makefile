# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../makefile.config

default: all

EXAMPLES = HelloMOAB GetEntities SetsNTags StructuredMeshSimple UniformRefinement DeleteEdges ReadPartFile ProjectToSphere
PAREXAMPLES = HelloParMOAB ReadWriteTest
EXOIIEXAMPLES = TestExodusII
ifeq ("$(MOAB_DEV)","yes")
ERROREXAMPLES = ErrorHandlingSimulation TestErrorHandling TestErrorHandlingPar
endif
ALLEXAMPLES = ${EXAMPLES}

ifeq ("$(MOAB_NETCDF_ENABLED)","yes")
ALLEXAMPLES += $(EXOIIEXAMPLES)
endif

ifeq ("$(MOAB_MPI_ENABLED)","yes")
parallel: ${PAREXAMPLES}
ALLEXAMPLES += ${PAREXAMPLES}
else
parallel:
endif

all: $(ALLEXAMPLES)

HelloMOAB: HelloMOAB.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-HelloMOAB: HelloMOAB
	${RUNSERIAL} ./HelloMOAB

ProjectToSphere: ProjectToSphere.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}
	
GetEntities: GetEntities.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

DeleteEdges: DeleteEdges.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-GetEntities: GetEntities
	${RUNSERIAL} ./GetEntities

run-DeleteEdges: DeleteEdges
	${RUNSERIAL} ./DeleteEdges

SetsNTags: SetsNTags.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-SetsNTags: SetsNTags
	${RUNSERIAL} ./SetsNTags

StructuredMeshSimple: StructuredMeshSimple.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-StructuredMeshSimple: StructuredMeshSimple
	${RUNSERIAL} ./StructuredMeshSimple -d 3 -n 5

ReadPartFile: ReadPartFile.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ReadPartFile: ReadPartFile
	${RUNSERIAL} ./ReadPartFile

HelloParMOAB: HelloParMOAB.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-HelloParMOAB: HelloParMOAB
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNPARALLEL} ./HelloParMOAB
endif

TestExodusII: TestExodusII.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-TestExodusII: TestExodusII
	${RUNSERIAL} ./TestExodusII

ReadWriteTest: ReadWriteTest.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ReadWriteTest: ReadWriteTest
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)-$(MOAB_NETCDF_ENABLED)","yes-yes-yes")
	${RUNPARALLEL} ./ReadWriteTest
endif

ErrorHandlingSimulation: ErrorHandlingSimulation.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestErrorHandling: TestErrorHandling.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

TestErrorHandlingPar: TestErrorHandlingPar.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

UniformRefinement: UniformRefinement.o
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

# -i is required; with no arguments ProjectToSphere tries to load "" and fails.
# A .vtk input keeps this independent of HDF5.
run-ProjectToSphere: ProjectToSphere
	${RUNSERIAL} ./ProjectToSphere -i ${MESH_DIR}/sphere_tris_5.vtk -o ProjectToSphere_out.vtk -R 1.0

run-UniformRefinement: UniformRefinement
ifeq ("$(MOAB_MPI_ENABLED)-$(MOAB_HDF5_ENABLED)","yes-yes")
	${RUNSERIAL}   ./UniformRefinement ${MESH_DIR}/64bricks_12ktet.h5m 2 2
	${RUNPARALLEL} ./UniformRefinement ${MESH_DIR}/64bricks_512hex_256part.h5m 2
else
	${RUNSERIAL} ./UniformRefinement ${MESH_DIR}/hex_2048.vtk 2
endif

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}


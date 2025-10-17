# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
MOAB_DIR := @CMAKE_INSTALL_PREFIX@

# MESH_DIR is the directory containing mesh files that come with MOAB source
MESH_DIR="@CMAKE_SOURCE_DIR@/MeshFiles/unittest"

MOAB_CMAKE="yes"

####### COMMON SETUP FOR ALL EXAMPLES ##########
ifneq ($(wildcard ${MOAB_DIR}/lib/moab.make),)

include ${MOAB_DIR}/lib/moab.make
RUNLD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${MOAB_DIR}/lib"

else

$(error Error cannot build examples without a valid MOAB_DIR (${MOAB_DIR}) build/installation path)

endif

default:

.SUFFIXES: .o .cpp .F90

VERBOSE=@
ifeq ($(V),1)
	VERBOSE=
endif

RUNSERIAL = LD_LIBRARY_PATH=${RUNLD_LIBRARY_PATH} 
ifeq ("$(MOAB_MPI_ENABLED)","yes")
NPROCS = @NP@
RUNPARALLEL = LD_LIBRARY_PATH=${RUNLD_LIBRARY_PATH} @MPIEXEC@ @MPIEXEC_NP@ @NP@
else
RUNPARALLEL = LD_LIBRARY_PATH=${RUNLD_LIBRARY_PATH} 
endif

.cpp.o:
	@echo "  [CXX]  $<"
	${VERBOSE}${MOAB_CXX} ${CXXFLAGS} ${MOAB_CXXFLAGS} ${MOAB_CPPFLAGS} ${MOAB_INCLUDES} -DMESH_DIR=\"${MESH_DIR}\" -c $<

.F90.o:
	@echo "   [FC]  $<"
  ${VERBOSE}${MOAB_FC} ${FCFLAGS} ${MOAB_CPPFLAGS} $(FC_DEFINE)MESH_DIR=\"${MESH_DIR}\" -c $<

info:
	@echo "Using installation MOAB_DIR = ${MOAB_DIR}"
	@echo "Using library paths during runs = ${RUNLD_LIBRARY_PATH}" 

clobber:
	@rm -rf *.o *.mod *.vtk


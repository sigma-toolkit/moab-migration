set(classes
  MoabReader)

include_directories(${PARAVIEW_INCLUDE_DIRS})
include_directories(.)

vtk_module_add_module(LMOAB
  CLASSES ${classes})

vtk_module_link(LMOAB
  PUBLIC
  MOAB)

paraview_add_server_manager_xmls(
   XMLS readers.xml)
   

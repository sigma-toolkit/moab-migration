#
# MOABTPLTargets.cmake
#
# moab_declare_tpl_target(<Namespaced::Name> <includes-var> <libraries-var>)
#
# MOAB's own Find modules predate imported targets: they set <PKG>_INCLUDES and
# <PKG>_LIBRARIES and nothing else.  Linking those raw values into libMOAB puts
# absolute paths straight into the exported link interface, which is what made
# an installed MOAB unusable on any machine whose third-party libraries sit
# somewhere else.  This wraps the two variables in an INTERFACE IMPORTED target
# so that libMOAB links a name, and MOABConfig.cmake re-creates the same name on
# the consumer side by re-running the Find module there.
#
# The target is GLOBAL so that subdirectories added after the find_package()
# call can link it, and is skipped if the caller already has one by that name -
# a consumer that found the same library through its own config package should
# keep its own definition.
#
function(moab_declare_tpl_target _target _includes_var _libraries_var)
  if(TARGET ${_target})
    return()
  endif()

  add_library(${_target} INTERFACE IMPORTED GLOBAL)
  if(${_includes_var})
    set_property(TARGET ${_target} PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES "${${_includes_var}}")
  endif()
  if(${_libraries_var})
    set_property(TARGET ${_target} PROPERTY
      INTERFACE_LINK_LIBRARIES "${${_libraries_var}}")
  endif()
endfunction()

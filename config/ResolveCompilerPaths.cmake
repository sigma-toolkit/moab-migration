#
# ResolveCompilerPaths - this module defines two macros
#
# ========================================================================================
# Macro: RESOLVE_LIBRARIES
#
# Description:
#   Parses a GNU Make-style link line, resolves all referenced libraries using CMake's
#   find_library(), and returns a de-duplicated, ordered list of full library paths.
#
#   It supports:
#     - "-l<name>": resolved using find_library()
#     - "-L<path>": adds to the search paths for subsequent -l entries
#     - Absolute paths to libraries: included directly if the file exists
#     - Skips duplicates
#     - Preserves order
#
# Parameters:
#   RESOLVED_LIBS_OUT (output)
#     The name of the output variable to store the resolved, ordered list of libraries.
#     Will contain full absolute paths to the libraries that were successfully resolved.
#
#   LINK_LINE (input)
#     A single string representing a GNU Make-style link line.
#     Example: "-L/foo/lib -lfoo -lbar /usr/lib/libbaz.a -lm"
#
# Usage Example:
#   RESOLVE_LIBRARIES("-L/opt/lib -lfoo -lbar -lm -lz" MY_RESOLVED_LIBS)
#   message(STATUS "Resolved libraries: ${MY_RESOLVED_LIBS}")
#
# Notes:
#   - Any unresolvable libraries will emit a warning.
#   - Non-library flags (like -pthread, -Wl, etc.) are ignored.
#   - Avoids using PARENT_SCOPE if macro is called from the top-level scope.
#
# ========================================================================================
#
# Macro: RESOLVE_INCLUDES
#
# Description:
#   Parses a compile line for `-I` include flags and extracts all valid include paths.
#   Only include paths that exist on the file system are returned. Duplicate paths are
#   automatically removed. Warnings are printed for non-existent paths.
#
# Parameters:
#   INCS (output)
#     The name of the variable to store the list of resolved include directories.
#
#   COMPILE_LINE (input)
#     A string representing a compile line (e.g., from pkg-config or compiler flags),
#     which may include -I flags (with or without quotes).
#     Example: "-I/usr/include -I\"/opt/libs/include\" -I./missing"
#
# Output:
#   - Sets the variable named by INCS with a list of valid, unique include paths.
#
# Usage Example:
#   RESOLVE_INCLUDES(MY_INCLUDE_DIRS "-I/usr/include -I\"/opt/special/include\" -I./fake")
#   message(STATUS "Resolved includes: ${MY_INCLUDE_DIRS}")
#
# ========================================================================================

macro(RESOLVE_LIBRARIES RESOLVED_LIBS_OUT LINK_LINE )
    # Clear output
    # Initialize the resolved libraries list to an empty list
    set(${RESOLVED_LIBS_OUT} "")
    
    # Create a list to keep track of libraries we have already resolved to avoid duplicates
    set(_seen_libs "")
    
    # Create a list to store library search paths defined by -L<path> flags
    set(_lib_search_paths "")

    # Check if LINK_LINE is a list (already resolved libraries)
    # Count the number of semicolons
    string(LENGTH "${LINK_LINE}" _len)
    string(REPLACE ";" "" _counted_var "${LINK_LINE}")
    string(LENGTH "${_counted_var}" _counted_var_len)
    math(EXPR _semicolon_count "${_len} - ${_counted_var_len}")
    if(_semicolon_count GREATER 0)
      separate_arguments(LINK_LINE_LIST UNIX_COMMAND "${LINK_LINE}")
      # If it's a list, we don't need to do any parsing. Just check each element.
      foreach(_lib IN LISTS LINK_LINE_LIST)
          # If it's not already in the resolved list, add it.
          list(FIND _seen_libs "${_lib}" _already_index)
          if(_already_index EQUAL -1)
            list(INSERT ${RESOLVED_LIBS_OUT} 0 "${_lib}")
            list(APPEND _seen_libs "${_lib}")
          endif()
      endforeach()

    else()
      # Tokenize the GNU-style link line into separate flags
      # Separate the link flags into individual components (e.g., -L/path, -lfoo, etc.)
      set(_link_flags "")
      separate_arguments(_link_flags UNIX_COMMAND "${LINK_LINE}")

          # Loop over each argument in the link flags list
      foreach(_flag IN LISTS _link_flags)
          # Handle -L<path> flag (library search directory)
          if(_flag MATCHES "^-L(.+)$")
              # Extract the directory path and add it to the search path list
              list(APPEND _lib_search_paths "${CMAKE_MATCH_1}")

          # Handle -l<library> flag (library to link)
          elseif(_flag MATCHES "^-l(.+)$")
              # Extract the library name (e.g., "foo" from -lfoo)
              set(_lib_name "${CMAKE_MATCH_1}")
              
              # Check if we have already processed this library to avoid duplicates
              list(FIND _seen_libs "${_lib_name}" _already_index)
              if(_already_index EQUAL -1)
                  # Library not seen before, we will try to resolve it
                  
                  unset(_lib_path CACHE)  # Clear any previous cache

                  # Platform-specific behavior for library suffixes
                  if(WIN32)
                      # On Windows, prefer to resolve ".lib" files for linking
                      find_library(_lib_path "${_lib_name}"
                          PATHS ${_lib_search_paths}
                          NO_DEFAULT_PATH
                          SUFFIXES "" ".lib"  # Only search for .lib files
                      )
                  else()
                      # On Linux/macOS, resolve shared libraries (.so/.dylib) and static libraries (.a)
                      find_library(_lib_path "${_lib_name}"
                          PATHS ${_lib_search_paths}
                          NO_DEFAULT_PATH
                      )
                  endif()

                  # Fallback to system-wide search if the library was not found in user-defined paths
                  if(NOT _lib_path)
                      find_library(_lib_path "${_lib_name}")
                  endif()

                  # If we found the library, add it to the resolved list
                  if(_lib_path)
                      list(APPEND ${RESOLVED_LIBS_OUT} "${_lib_path}")
                      list(APPEND _seen_libs "${_lib_name}")
                  else()
                      # Print a warning if the library could not be found
                      message(WARNING "Library '${_lib_name}' not found.")
                  endif()
              endif()

          # Handle absolute paths to libraries (e.g., /path/to/libfoo.a)
          elseif(IS_ABSOLUTE "${_flag}" AND EXISTS "${_flag}")
              # Add the library path directly to the resolved list if it exists
              list(FIND _seen_libs "${_flag}" _already_index)
              if(_already_index EQUAL -1)
                  # Optional: Warn if a DLL is passed as a library path
                  if(WIN32 AND _flag MATCHES "\\.dll$")
                      message(WARNING "Ignoring DLL file for linking: ${_flag}")
                  else()
                      list(APPEND ${RESOLVED_LIBS_OUT} "${_flag}")
                  endif()
                  list(APPEND _seen_libs "${_flag}")
              endif()

          # If the flag is not recognized, print a status message (for debugging)
          else()
              message(STATUS "Ignoring unrecognized flag: ${_flag}")
          endif()
      endforeach()
    endif()

    # Return the resolved libraries list to the caller 
    set(${RESOLVED_LIBS_OUT} "${${RESOLVED_LIBS_OUT}}")
endmacro()



macro (RESOLVE_INCLUDES_OLD INCS COMPILE_LINE)
  string (REGEX MATCHALL "-I([^\" ]+|\"[^\"]+\")" _all_tokens "${COMPILE_LINE}")
  set (_incs_found)
  foreach (token ${_all_tokens})
    string (REGEX REPLACE "^-I" "" token ${token})
    string (REGEX REPLACE "//" "/" token ${token})
    if (EXISTS ${token})
      list (APPEND _incs_found ${token})
    else (EXISTS ${token})
      message (STATUS "Include directory ${token} does not exist")
    endif (EXISTS ${token})
  endforeach (token)
  if (_incs_found)
    list (REMOVE_DUPLICATES _incs_found)
  endif(_incs_found)
  set (${INCS} "${_incs_found}")
endmacro (RESOLVE_INCLUDES_OLD)

macro(RESOLVE_INCLUDES INCS COMPILE_LINE)
    string(REGEX MATCHALL "-I([^\" ]+|\"[^\"]+\")" _all_tokens "${COMPILE_LINE}")
    set(_incs_found)

    foreach(token ${_all_tokens})
        string(REGEX REPLACE "^-I" "" token "${token}")
        string(REGEX REPLACE "//" "/" token "${token}")

        # Windows: Normalize backslashes
        if(WIN32)
            string(REPLACE "\\" "/" token "${token}")
        endif()

        # Remove surrounding quotes if any
        string(REGEX REPLACE "^\"(.*)\"$" "\\1" token "${token}")

        if(EXISTS "${token}")
            list(APPEND _incs_found "${token}")
        else()
            message(STATUS "Include directory '${token}' does not exist")
        endif()
    endforeach()

    if(_incs_found)
        list(REMOVE_DUPLICATES _incs_found)
    endif()

    set(${INCS} "${_incs_found}")
endmacro()

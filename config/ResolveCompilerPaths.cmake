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

# ----------------------------------------------------------------------------
# Linker-token exception table (shared by RESOLVE_LIBRARIES and
# FILTER_LINK_LIBRARIES).  Items matching any of these regular expressions
# are kept verbatim instead of being treated as unrecognized flags.  This is
# required for platform-specific linker syntax that find_library cannot
# resolve (e.g. Apple "-framework Accelerate", GNU "-Wl,...", "-Xlinker ...").
# Extend these lists when a new linker construct is being incorrectly dropped.
# ----------------------------------------------------------------------------
set(_RCP_KEEP_PATTERNS
    "^-framework"     # macOS frameworks: "-framework Accelerate" (combined or split)
    "^-Wl,"           # GNU linker pass-through flags
    "^-Xlinker"       # Alternative linker pass-through
    "^-pthread$"      # POSIX threads
    )

# Subset of keep-patterns whose match consumes the next list element as its
# argument (e.g. "-framework Accelerate" given as two tokens).  Both tokens
# must be preserved together and in order.
set(_RCP_FLAGS_WITH_ARG "^-framework$|^-Xlinker$")

# Internal helper: returns TRUE in OUT_VAR if TOKEN matches any keep-pattern,
# and TRUE in OUT_CONSUMES if it also consumes the following token.
function(_RCP_TOKEN_IS_KEPT TOKEN OUT_VAR OUT_CONSUMES)
    set(${OUT_VAR} FALSE PARENT_SCOPE)
    set(${OUT_CONSUMES} FALSE PARENT_SCOPE)
    foreach(_pat IN LISTS _RCP_KEEP_PATTERNS)
        if(TOKEN MATCHES "${_pat}")
            set(${OUT_VAR} TRUE PARENT_SCOPE)
            if(TOKEN MATCHES "${_RCP_FLAGS_WITH_ARG}")
                set(${OUT_CONSUMES} TRUE PARENT_SCOPE)
            endif()
            return()
        endif()
    endforeach()
endfunction()


# ============================================================================
# Macro: FILTER_LINK_LIBRARIES(<listvar>)
#
# Validates a CMake list of link-line items in-place, dropping anything that
# is neither (a) an absolute path to an existing file, (b) a -l flag, (c) a
# plain library name (e.g. "stdc++"), nor (d) a token covered by the keep-
# pattern table above.  Items that consume the next argument (e.g.
# "-framework") keep their following token verbatim and in order.
#
# Use this on link-library lists assembled from heterogeneous sources (find
# modules, hand-set BLAS/LAPACK, etc.) before handing them to
# target_link_libraries().  Items dropped are reported via message(STATUS).
# ============================================================================
macro(FILTER_LINK_LIBRARIES _var)
    set(_flb_out "")
    set(_flb_pending "")     # holds a flag like "-framework" awaiting its arg
    foreach(_flb_item IN LISTS ${_var})
        # Pair-completion: combine a pending consuming flag with its argument
        # into a single space-joined string (e.g. "-framework Accelerate"),
        # so that target_link_libraries() does not mistake the bare argument
        # for a library name and emit "-l<arg>".
        if(NOT _flb_pending STREQUAL "")
            list(APPEND _flb_out "${_flb_pending} ${_flb_item}")
            set(_flb_pending "")
            continue()
        endif()

        _RCP_TOKEN_IS_KEPT("${_flb_item}" _flb_keep _flb_consumes)
        if(_flb_keep)
            if(_flb_consumes)
                # Defer until we see the following token to join them.
                set(_flb_pending "${_flb_item}")
            else()
                list(APPEND _flb_out "${_flb_item}")
            endif()
            continue()
        endif()

        if(IS_ABSOLUTE "${_flb_item}")
            if(EXISTS "${_flb_item}" AND NOT IS_DIRECTORY "${_flb_item}")
                list(APPEND _flb_out "${_flb_item}")
            else()
                message(STATUS "Dropping non-file path from link libraries: ${_flb_item}")
            endif()
        elseif(_flb_item MATCHES "^-l")
            list(APPEND _flb_out "${_flb_item}")
        elseif(_flb_item MATCHES "^[A-Za-z]")
            # Plain library name (e.g. "stdc++", "m", "z") — pass through
            list(APPEND _flb_out "${_flb_item}")
        else()
            message(STATUS "Dropping unrecognized item from link libraries: ${_flb_item}")
        endif()
    endforeach()

    # A trailing consuming flag with no following token is malformed; warn.
    if(NOT _flb_pending STREQUAL "")
        message(WARNING "Dangling consuming flag '${_flb_pending}' with no argument in ${_var}")
    endif()

    set(${_var} ${_flb_out})
    unset(_flb_out)
    unset(_flb_item)
    unset(_flb_keep)
    unset(_flb_consumes)
    unset(_flb_pending)
endmacro()


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
      # Pre-tokenized input.  Preserve original order (use APPEND, not
      # INSERT 0) and honor the keep-pattern table.  Consuming flags such
      # as "-framework Accelerate" are emitted as a single space-joined
      # element so downstream target_link_libraries() calls do not parse
      # the bare argument as a library name.
      set(_rl_pending "")
      foreach(_lib IN LISTS LINK_LINE_LIST)
          if(NOT _rl_pending STREQUAL "")
              set(_rl_combined "${_rl_pending} ${_lib}")
              list(FIND _seen_libs "${_rl_combined}" _already_index)
              if(_already_index EQUAL -1)
                  list(APPEND ${RESOLVED_LIBS_OUT} "${_rl_combined}")
                  list(APPEND _seen_libs "${_rl_combined}")
              endif()
              set(_rl_pending "")
              continue()
          endif()

          _RCP_TOKEN_IS_KEPT("${_lib}" _rl_keep _rl_consumes)
          if(_rl_keep)
              if(_rl_consumes)
                  set(_rl_pending "${_lib}")
              else()
                  list(FIND _seen_libs "${_lib}" _already_index)
                  if(_already_index EQUAL -1)
                      list(APPEND ${RESOLVED_LIBS_OUT} "${_lib}")
                      list(APPEND _seen_libs "${_lib}")
                  endif()
              endif()
              continue()
          endif()

          # If it's not already in the resolved list, add it.
          list(FIND _seen_libs "${_lib}" _already_index)
          if(_already_index EQUAL -1)
              list(APPEND ${RESOLVED_LIBS_OUT} "${_lib}")
              list(APPEND _seen_libs "${_lib}")
          endif()
      endforeach()
      if(NOT _rl_pending STREQUAL "")
          message(WARNING "Dangling consuming flag '${_rl_pending}' with no argument in input: ${LINK_LINE}")
      endif()

    else()
      # Tokenize the GNU-style link line into separate flags
      # Separate the link flags into individual components (e.g., -L/path, -lfoo, etc.)
      set(_link_flags "")
      separate_arguments(_link_flags UNIX_COMMAND "${LINK_LINE}")

          # Loop over each argument in the link flags list.  Consuming
          # flags (e.g. "-framework Accelerate") are deferred so they can
          # be emitted as a single space-joined element — keeping the pair
          # atomic for downstream target_link_libraries() consumers.
      set(_rl_pending "")
      foreach(_flag IN LISTS _link_flags)
          if(NOT _rl_pending STREQUAL "")
              set(_rl_combined "${_rl_pending} ${_flag}")
              list(FIND _seen_libs "${_rl_combined}" _already_index)
              if(_already_index EQUAL -1)
                  list(APPEND ${RESOLVED_LIBS_OUT} "${_rl_combined}")
                  list(APPEND _seen_libs "${_rl_combined}")
              endif()
              set(_rl_pending "")
              continue()
          endif()

          # Pass through platform-specific linker syntax (frameworks, -Wl,
          # -Xlinker, -pthread, ...) covered by the shared keep-pattern table.
          _RCP_TOKEN_IS_KEPT("${_flag}" _rl_keep _rl_consumes)
          if(_rl_keep)
              if(_rl_consumes)
                  set(_rl_pending "${_flag}")
              else()
                  list(FIND _seen_libs "${_flag}" _already_index)
                  if(_already_index EQUAL -1)
                      list(APPEND ${RESOLVED_LIBS_OUT} "${_flag}")
                      list(APPEND _seen_libs "${_flag}")
                  endif()
              endif()
              continue()
          endif()

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
      if(NOT _rl_pending STREQUAL "")
          message(WARNING "Dangling consuming flag '${_rl_pending}' with no argument in input: ${LINK_LINE}")
      endif()
    endif()

    # Return the resolved libraries list to the caller
    set(${RESOLVED_LIBS_OUT} "${${RESOLVED_LIBS_OUT}}")

    # Tidy up macro-local helper variables to avoid leaking into the caller.
    unset(_rl_pending)
    unset(_rl_combined)
    unset(_rl_keep)
    unset(_rl_consumes)
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

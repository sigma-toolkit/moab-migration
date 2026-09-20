function(_moab_add_dist_targets)
  if(NOT UNIX)
    message(WARNING "The dist targets currently support only UNIX-like platforms.")
    return()
  endif()

  find_program(GIT_EXECUTABLE git REQUIRED)
  find_program(GPG_EXECUTABLE gpg)
  find_program(PERL_EXECUTABLE perl REQUIRED)
  option(MOAB_DIST_SIGN "Sign distribution archives with GPG" OFF)
  if(MOAB_DIST_SIGN AND NOT GPG_EXECUTABLE)
    message(FATAL_ERROR "MOAB_DIST_SIGN requires gpg")
  endif()
  set(_dist_sign_args)
  set(_dist_sign_comment "")
  if(MOAB_DIST_SIGN)
    list(APPEND _dist_sign_args -DGPG_EXECUTABLE=${GPG_EXECUTABLE})
    set(_dist_sign_comment " and signature")
  endif()

  set(MOAB_DIST_NAME "moab-${CMAKE_PROJECT_VERSION}" PARENT_SCOPE)
  set(MOAB_DIST_DIR "${CMAKE_BINARY_DIR}/moab-${CMAKE_PROJECT_VERSION}" PARENT_SCOPE)
  set(MOAB_DIST_ARTIFACT_FILE "${CMAKE_BINARY_DIR}/dist-tarball-path.txt" PARENT_SCOPE)
  set(MOAB_DIST_BASENAME_FILE "${CMAKE_BINARY_DIR}/dist-archive-basename.txt" PARENT_SCOPE)
  set(_dist_script_args
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
    -DBINARY_DIR=${CMAKE_BINARY_DIR}
    -DPROJECT_VERSION=${CMAKE_PROJECT_VERSION}
    -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
    -DPERL_EXECUTABLE=${PERL_EXECUTABLE}
    -P "${CMAKE_SOURCE_DIR}/config/dist.cmake")

  add_custom_target(distdir
    COMMAND "${CMAKE_COMMAND}" -DMODE=distdir ${_dist_script_args}
    COMMENT "Generating distribution directory from tracked working-tree files..."
    VERBATIM)

  foreach(format IN ITEMS gz bz2 xz)
    if(format STREQUAL "gz")
      set(target dist_targz)
    elseif(format STREQUAL "bz2")
      set(target dist_tarbz2)
    else()
      set(target dist_tarxz)
    endif()
    add_custom_target(${target}
      COMMAND "${CMAKE_COMMAND}" -DMODE=archive -DARCHIVE_FORMAT=${format}
        ${_dist_sign_args} ${_dist_script_args}
      COMMENT "Generating provenance-named tar.${format} archive${_dist_sign_comment}..."
      VERBATIM)
    add_dependencies(${target} distdir)
  endforeach()

  add_custom_target(dist DEPENDS dist_targz)
  add_custom_target(distorig
    COMMAND "${CMAKE_COMMAND}" -DMODE=copy_orig ${_dist_script_args}
    VERBATIM)
  add_dependencies(distorig dist)
  add_custom_target(distclean
    COMMAND "${CMAKE_COMMAND}" -DMODE=clean ${_dist_script_args}
    COMMENT "Cleaning distribution sources and archives..."
    VERBATIM)
endfunction()

if(MODE)
  if(NOT SOURCE_DIR OR NOT BINARY_DIR OR NOT PROJECT_VERSION OR NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "The distribution script is missing a required argument")
  endif()

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE git_result
    OUTPUT_VARIABLE git_hash
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT git_result EQUAL 0)
    message(FATAL_ERROR "Unable to determine the HEAD commit hash")
  endif()
  string(SUBSTRING "${git_hash}" 0 12 short_hash)

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" diff --quiet HEAD --
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE dirty_result)
  if(dirty_result EQUAL 0)
    set(dirty_suffix "")
    set(dirty_state "clean")
  elseif(dirty_result EQUAL 1)
    set(dirty_suffix "-dirty")
    set(dirty_state "dirty (tracked working-tree content differs from HEAD)")
  else()
    message(FATAL_ERROR "Unable to determine the tracked working-tree state")
  endif()

  set(dist_name "moab-${PROJECT_VERSION}")
  set(dist_dir "${BINARY_DIR}/${dist_name}")
  set(archive_base "${dist_name}-${short_hash}${dirty_suffix}")
  set(artifact_file "${BINARY_DIR}/dist-tarball-path.txt")
  set(basename_file "${BINARY_DIR}/dist-archive-basename.txt")
endif()

if(MODE STREQUAL "distdir")
  if(NOT PERL_EXECUTABLE)
    message(FATAL_ERROR "The distdir script requires Perl")
  endif()

  set(base_archive "${BINARY_DIR}/dist-head.tar")
  set(working_tree_patch "${BINARY_DIR}/dist-working-tree.patch")
  file(REMOVE_RECURSE "${dist_dir}")
  file(REMOVE "${base_archive}" "${working_tree_patch}")
  file(MAKE_DIRECTORY "${dist_dir}")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" archive --format=tar --output "${base_archive}" HEAD
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE archive_result)
  if(NOT archive_result EQUAL 0)
    message(FATAL_ERROR "Unable to export the HEAD source tree")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E tar xf "${base_archive}"
    WORKING_DIRECTORY "${dist_dir}"
    RESULT_VARIABLE extract_result)
  if(NOT extract_result EQUAL 0)
    message(FATAL_ERROR "Unable to extract the HEAD source tree")
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" diff --binary --full-index HEAD --
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE diff_result
    OUTPUT_FILE "${working_tree_patch}")
  if(NOT diff_result EQUAL 0)
    message(FATAL_ERROR "Unable to capture tracked working-tree changes")
  endif()
  file(SIZE "${working_tree_patch}" patch_size)
  if(patch_size GREATER 0)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" apply --unsafe-paths --directory "${dist_dir}" "${working_tree_patch}"
      WORKING_DIRECTORY "${SOURCE_DIR}"
      RESULT_VARIABLE apply_result)
    if(NOT apply_result EQUAL 0)
      message(FATAL_ERROR "Unable to apply tracked working-tree changes to the distribution")
    endif()
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE final_hash_result
    OUTPUT_VARIABLE final_git_hash
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" diff --binary --full-index HEAD --
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE final_diff_result
    OUTPUT_VARIABLE final_patch)
  file(READ "${working_tree_patch}" initial_patch)
  if(NOT final_hash_result EQUAL 0 OR NOT final_diff_result EQUAL 0 OR
     NOT "${git_hash}" STREQUAL "${final_git_hash}" OR NOT "${initial_patch}" STREQUAL "${final_patch}")
    message(FATAL_ERROR "Tracked sources changed while the distribution snapshot was being created")
  endif()
  file(REMOVE "${base_archive}" "${working_tree_patch}")

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" log -1 --format=%B HEAD
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE message_result
    OUTPUT_VARIABLE commit_message
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT message_result EQUAL 0)
    message(FATAL_ERROR "Unable to determine the HEAD commit message")
  endif()
  file(WRITE "${dist_dir}/GITHASH"
    "HEAD: ${git_hash}\nState: ${dirty_state}\n\nCommit message:\n${commit_message}\n")
  file(WRITE "${basename_file}" "${archive_base}\n")
  file(WRITE "${dist_dir}/.version" "${PROJECT_VERSION}\n")
  execute_process(
    COMMAND "${PERL_EXECUTABLE}" "${SOURCE_DIR}/config/gitlog-to-changelog"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE changelog_result
    OUTPUT_FILE "${dist_dir}/ChangeLog")
  if(NOT changelog_result EQUAL 0)
    message(FATAL_ERROR "Unable to generate ChangeLog")
  endif()
elseif(MODE STREQUAL "archive")
  if(NOT EXISTS "${dist_dir}" OR NOT EXISTS "${basename_file}")
    message(FATAL_ERROR "Run the distdir target before creating an archive")
  endif()
  file(READ "${basename_file}" archive_base)
  string(STRIP "${archive_base}" archive_base)
  if(ARCHIVE_FORMAT STREQUAL "gz")
    set(extension "tar.gz")
    set(tar_flag czf)
  elseif(ARCHIVE_FORMAT STREQUAL "bz2")
    set(extension "tar.bz2")
    set(tar_flag cjf)
  elseif(ARCHIVE_FORMAT STREQUAL "xz")
    set(extension "tar.xz")
    set(tar_flag cJf)
  else()
    message(FATAL_ERROR "Unsupported archive format: ${ARCHIVE_FORMAT}")
  endif()

  set(archive "${BINARY_DIR}/${archive_base}.${extension}")
  file(REMOVE "${archive}" "${archive}.sig")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E tar "${tar_flag}" "${archive}" "${dist_name}"
    WORKING_DIRECTORY "${BINARY_DIR}"
    RESULT_VARIABLE archive_result)
  if(NOT archive_result EQUAL 0)
    message(FATAL_ERROR "Unable to create ${archive}")
  endif()
  if(GPG_EXECUTABLE)
    execute_process(
      COMMAND "${GPG_EXECUTABLE}" --detach-sign --armor -o "${archive}.sig" "${archive}"
      RESULT_VARIABLE signature_result)
    if(NOT signature_result EQUAL 0)
      message(FATAL_ERROR "Unable to sign ${archive}")
    endif()
  endif()
  if(ARCHIVE_FORMAT STREQUAL "gz")
    file(WRITE "${artifact_file}" "${archive}\n")
  endif()
  message(STATUS "Created ${archive}")
elseif(MODE STREQUAL "copy_orig")
  if(NOT EXISTS "${artifact_file}")
    message(FATAL_ERROR "Run the dist target before distorig")
  endif()
  file(READ "${artifact_file}" archive)
  string(STRIP "${archive}" archive)
  file(COPY_FILE "${archive}" "${BINARY_DIR}/${dist_name}.orig.tar.gz")
elseif(MODE STREQUAL "clean")
  file(GLOB archives
    "${BINARY_DIR}/${dist_name}-*.tar.gz"
    "${BINARY_DIR}/${dist_name}-*.tar.gz.sig"
    "${BINARY_DIR}/${dist_name}-*.tar.bz2"
    "${BINARY_DIR}/${dist_name}-*.tar.bz2.sig"
    "${BINARY_DIR}/${dist_name}-*.tar.xz"
    "${BINARY_DIR}/${dist_name}-*.tar.xz.sig")
  file(REMOVE_RECURSE "${dist_dir}")
  file(REMOVE "${artifact_file}" "${basename_file}" "${BINARY_DIR}/${dist_name}.orig.tar.gz" ${archives})
endif()

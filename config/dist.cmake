MACRO(_SETUP_PROJECT_DIST)
  IF(UNIX)
    FIND_PROGRAM(TAR tar)
    FIND_PROGRAM(GPG gpg)

    IF(APPLE)
      SET(IS_MACOS TRUE)
    ELSE()
      SET(IS_MACOS FALSE)
    ENDIF()

    SET(DIST_TMP_FILE ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz)

    ADD_CUSTOM_TARGET(distdir
      COMMAND
      ${CMAKE_COMMAND} -E remove -f ${DIST_TMP_FILE}
      #&& ${CMAKE_SOURCE_DIR}/config/git-archive-all.sh --prefix=${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION} ${DIST_TMP_FILE}
      #&& git archive --prefix=${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}/ ${DIST_TMP_FILE}
      && git archive -o ${DIST_TMP_FILE} --format=tar.gz --prefix=${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}/ HEAD "\":(exclude)deprecated/\"" "\":(exclude).circleci/\"" "\":(exclude)debian/\""
      #git archive -v -o eb-bundle.zip --format=zip HEAD . ":(exclude)data/local.js"
      && cd ${CMAKE_BINARY_DIR}
      && ${TAR} -xf ${DIST_TMP_FILE}
      && echo "${CMAKE_PROJECT_VERSION}" > ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}/.version
      && ${CMAKE_SOURCE_DIR}/config/gitlog-to-changelog > ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}/ChangeLog
      && ${CMAKE_COMMAND} -E remove -f ${DIST_TMP_FILE}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      COMMENT "Generating dist directory..."
    )

    ADD_CUSTOM_TARGET(dist_targz
      COMMAND
      ${TAR} -czf ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}
      && ${GPG} --detach-sign --armor -o ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz.sig ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generating tar.gz tarball and its signature..."
    )

    ADD_CUSTOM_TARGET(dist_tarbz2
      COMMAND
      ${TAR} -cjf ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.bz2 ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}
      && ${GPG} --detach-sign --armor -o ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.bz2.sig ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.bz2
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generating tar.bz2 tarball and its signature..."
    )

    ADD_CUSTOM_TARGET(dist_tarxz
      COMMAND
      ${TAR} -cJf ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.xz ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}
      && ${GPG} --detach-sign --armor -o ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.xz.sig ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.xz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generating tar.xz tarball and its signature..."
    )

    ADD_CUSTOM_TARGET(dist DEPENDS dist_targz)

    ADD_CUSTOM_TARGET(distclean
      COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Cleaning dist sources..."
    )

    ADD_CUSTOM_TARGET(distorig
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.orig.tar.gz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generating orig tarball..."
    )

    ADD_DEPENDENCIES(dist_targz distdir)
    ADD_DEPENDENCIES(dist_tarbz2 distdir)
    ADD_DEPENDENCIES(dist_tarxz distdir)
    ADD_DEPENDENCIES(distorig dist)
  ELSE()
    MESSAGE(WARNING "The _SETUP_PROJECT_DIST macro currently supports only UNIX-like platforms.")
  ENDIF()
ENDMACRO(_SETUP_PROJECT_DIST)

# Macro: TEST_BIG_ENDIAN(VARIABLE)
# Purpose:
#   Detects the endianness of the system (big endian or little endian) and
#   stores the result in the provided variable.
#
# Arguments:
#   VARIABLE - The variable to store the result into:
#              0 = little endian
#              1 = big endian
#
# Notes:
#   - Uses a small test program to write a known pattern to a binary file,
#     and inspects the result to determine byte order.
#   - Supports CMake 3.x and 4.x.

macro(TEST_BIG_ENDIAN VARIABLE)
  if(HAVE_${VARIABLE} MATCHES ^HAVE_${VARIABLE}$)
    message(STATUS "Checking if the system is big endian")

    # Include the check_type_size macro to find a suitable 16-bit type
    include(CheckTypeSize)

    # Find a usable 16-bit unsigned integer type
    check_type_size("unsigned short" CMAKE_SIZEOF_UNSIGNED_SHORT)
    if(CMAKE_SIZEOF_UNSIGNED_SHORT EQUAL 2)
      set(CMAKE_16BIT_TYPE "unsigned short")
    else()
      check_type_size("unsigned int" CMAKE_SIZEOF_UNSIGNED_INT)
      if(CMAKE_SIZEOF_UNSIGNED_INT EQUAL 2)
        set(CMAKE_16BIT_TYPE "unsigned int")
      else()
        check_type_size("unsigned long" CMAKE_SIZEOF_UNSIGNED_LONG)
        if(CMAKE_SIZEOF_UNSIGNED_LONG EQUAL 2)
          set(CMAKE_16BIT_TYPE "unsigned long")
        else()
          message(FATAL_ERROR "No suitable 16-bit type found")
        endif()
      endif()
    endif()

    # Write a small C program to test endianness
    set(TEST_ENDIANESS_OUTPUT_C "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/TestEndianess.c")
    file(WRITE "${TEST_ENDIANESS_OUTPUT_C}" "
#include <stdio.h>
typedef ${CMAKE_16BIT_TYPE} cmakeint16;

int main() {
  cmakeint16 val = 0x0102;
  FILE *f = fopen(\"TestEndianess.bin\", \"wb\");
  if (!f) return 1;
  if (*((unsigned char*)&val) == 0x02)
    fputs(\"THIS IS LITTLE ENDIAN\\n\", f);
  else
    fputs(\"THIS IS BIG ENDIAN\\n\", f);
  fclose(f);
  return 0;
}
")

    # Keep the source content to log later in case of failure
    file(READ "${TEST_ENDIANESS_OUTPUT_C}" TEST_ENDIANESS_FILE_CONTENT)

    # Try to compile and run the test program
    try_compile(HAVE_${VARIABLE}
      "${CMAKE_BINARY_DIR}"
      "${TEST_ENDIANESS_OUTPUT_C}"
      OUTPUT_VARIABLE OUTPUT
      COPY_FILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestEndianess.bin")

    if(HAVE_${VARIABLE})
      # Read output from binary file to detect endianness
      file(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestEndianess.bin"
           CMAKE_TEST_ENDIANESS_STRINGS_LE LIMIT_COUNT 1 REGEX "THIS IS LITTLE ENDIAN")

      file(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestEndianess.bin"
           CMAKE_TEST_ENDIANESS_STRINGS_BE LIMIT_COUNT 1 REGEX "THIS IS BIG ENDIAN")

      # Handle edge case for universal binaries (e.g., on macOS)
      if(CMAKE_TEST_ENDIANESS_STRINGS_BE AND CMAKE_TEST_ENDIANESS_STRINGS_LE)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "powerpc")
          set(CMAKE_TEST_ENDIANESS_STRINGS_BE TRUE)
          set(CMAKE_TEST_ENDIANESS_STRINGS_LE FALSE)
        else()
          set(CMAKE_TEST_ENDIANESS_STRINGS_BE FALSE)
          set(CMAKE_TEST_ENDIANESS_STRINGS_LE TRUE)
        endif()
        message(STATUS "TEST_BIG_ENDIAN found conflicting results on a universal binary system.")
      endif()

      # Set result
      if(CMAKE_TEST_ENDIANESS_STRINGS_LE)
        set(${VARIABLE} 0 CACHE INTERNAL "Result of TEST_BIG_ENDIAN" FORCE)
        message(STATUS "System is little endian")
      endif()

      if(CMAKE_TEST_ENDIANESS_STRINGS_BE)
        set(${VARIABLE} 1 CACHE INTERNAL "Result of TEST_BIG_ENDIAN" FORCE)
        message(STATUS "System is big endian")
      endif()

      if(NOT CMAKE_TEST_ENDIANESS_STRINGS_BE AND NOT CMAKE_TEST_ENDIANESS_STRINGS_LE)
        message(SEND_ERROR "TEST_BIG_ENDIAN found no match in compiled binary!")
      endif()

      # Append full output and source to error log for diagnostics
      file(APPEND "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log"
           "Determining endianness passed:\n${OUTPUT}\nTestEndianess.c:\n${TEST_ENDIANESS_FILE_CONTENT}\n\n")
    else()
      # Compilation failed
      message(STATUS "Endian test compilation failed")
      file(APPEND "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log"
           "Determining endianness failed:\n${OUTPUT}\nTestEndianess.c:\n${TEST_ENDIANESS_FILE_CONTENT}\n\n")
      set(${VARIABLE})
    endif()
  endif()
endmacro()

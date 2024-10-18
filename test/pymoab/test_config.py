import subprocess
import pytest

def test_cmake_build(tmp_path):
    """
    Test if CMake can find the MOAB package correctly and build a simple test library.
    """
    # Define the path to the temporary build directory.
    build_dir = tmp_path / "build"
    build_dir.mkdir()
    
    # Create a minimal CPP file in the temporary directory.
    cpp_content = """
    #include <moab/Core.hpp>
    int main() {
        return 0; // Basic test to ensure MOAB can be included.
    }
    """
    cpp_file = tmp_path / "test.cpp"
    cpp_file.write_text(cpp_content)

    # Create a minimal CMakeLists.txt file in the temporary directory.
    cmake_content = """
    cmake_minimum_required(VERSION 3.15)
    project(TestMOAB)
    find_package(MOAB REQUIRED)
    include_directories(${MOAB_INCLUDE_DIRS})
    add_executable(test_moab test.cpp)
    target_link_libraries(test_moab PUBLIC MOAB)
    """
    cmake_file = tmp_path / "CMakeLists.txt"
    cmake_file.write_text(cmake_content)

    # Run CMake in the temporary build directory.
    try:
        subprocess.run(
            ["cmake", str(tmp_path), "-DMOAB_DIR=$MOAB_DIR"],
            cwd=build_dir,
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        pytest.fail(f"CMake configuration failed with error: {e.stderr}")

    # Build the project using CMake.
    try:
        subprocess.run(
            ["cmake", "--build", "."],
            cwd=build_dir,
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        pytest.fail(f"CMake build failed with error: {e.stderr}")

    # If the build is successful, the test will pass.
    test_executable = build_dir / "test_moab"
    assert test_executable.exists(), "The test executable was not created."

import os
from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake, cmake_layout

class Moab(ConanFile):
    version = "5.5.1"

    name = "moab"
    license = "LGPL-3.0/BSD"
    author = "moab-dev@anl.gov"
    url = "https://bitbucket.org/fathomteam/moab/src/master/"
    description = "MOAB: Mesh-Oriented datABase"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False], "plugin": [True, False], "hdf5": [True, False], "parallel": [True, False]}
    default_options = {"shared": True, "fPIC": True, "plugin": False, "hdf5": True, "parallel": False}
#    generators = "CMakeToolchain", "CMakeDeps"

    exports_sources = "*", "!build/*"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def requirements(self):
        self.requires("eigen/3.3.9")
        if self.options.hdf5:
            self.requires("hdf5/1.12.0", options={"shared": True})
        if self.options.parallel:
            self.requires("openmpi/4.1.0")
        self.requires("zlib/1.2.11")
        self.requires("libcurl/8.10.1")
        if self.options.plugin:
            self.requires("paraview/5.12.1@d3d/stable")

    def layout(self):
        # We make the assumption that if the compiler is msvc the
        # CMake generator is multi-config
        multi = True if self.settings.get_safe("compiler") == "msvc" else False
        if multi:
            self.folders.generators = os.path.join("build", "generators")
            self.folders.build = "build"
        else:
            if self.settings.build_type == "Debug":
                self.folders.generators = os.path.join("dbg", "cmake")
                self.folders.build = "dbg"
            else:
                self.folders.generators = os.path.join("opt", "cmake")
                self.folders.build = "opt"

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["ENABLE_FORTRAN"] = False
        if self.options.hdf5:
            tc.variables["ENABLE_HDF5"] = True
        if self.options.parallel:
            tc.variables["ENABLE_MPI"] = True
        tc.variables["ENABLE_ZLIB"] = True
        tc.variables["ENABLE_EIGEN3"] = True
        tc.variables["ENABLE_BLASLAPACK"] = False
        tc.variables["ENABLE_IMESH"] = False
        tc.variables["ENABLE_VTK"] = self.options.plugin
        tc.variables["MOAB_BUILD_VTKMOABREADER"] = self.options.plugin
        tc.generate()
        tc2 = CMakeDeps(self)
        tc2.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.env_info.path.append(os.path.join(self.package_folder, "bin"))
        self.cpp_info.libs = ["MOAB"]

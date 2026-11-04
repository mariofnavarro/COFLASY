from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout

class CoflasyCConan(ConanFile):
    name = "coflasy-c"
    version = "0.1"
    settings = "os", "compiler", "build_type", "arch"

    requires = (
        "gsl/2.7.1"
    )
    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeDeps(self).generate()
        CMakeToolchain(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
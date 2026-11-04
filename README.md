# COFLASY-C: COsmological evolution of FLAvour ASYmmetries - C++

The C++ code "COFLASY-C" has been developed by Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner in order to study the evolution of primordial lepton flavour asymmetries from the early Universe to BBN, following the approach of arXiv:2502.14960. If you use this code for scientific publications, please cite these papers: 

"Lepton Flavour Asymmetries: from the early Universe to BBN",  Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner,
[JHEP 06 (2025), 137] ( https://link.springer.com/article/10.1007/JHEP06(2025)137 ), [arXiv:2502.14960] ( https://arxiv.org/abs/2502.14960 ).

"A Limit on the Total Lepton Number in the Universe from BBN and the CMB", Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner,
[JCAP 02 (2026) 017] ( https://iopscience.iop.org/article/10.1088/1475-7516/2026/02/017 ), [arXiv:2510.02438] ( https://arxiv.org/abs/2510.02438 ).

## Requirements

| Package | Description |
| ------ | ----------- |
| [CMake](https://cmake.org)   | generate the Makefile |
| [GSL](https://www.gnu.org/software/gsl/) | numerical library for C++ |
| [Conan](https://conan.io) (optional) | C++ package manager; fetches GSL and builds the code |


## Installation

### Option 1: CMake

With GSL installed on your system, the Makefile can be generated within the main folder by executing the following:

    mkdir build
    cd build
    cmake ..
    make

This should create the executable coflasy-c.exe inside the build directory.

### Option 2: Conan (recommended)

With [Conan 2](https://conan.io) installed, the code can be built from the main folder with a single command. Conan fetches the GSL dependency automatically, so it does not need to be installed by hand:

    conan build .

If this is your first time using Conan, create a default profile first with `conan profile detect`. If no prebuilt GSL binary is available for your system, add `--build=missing` so that Conan builds it from source:

    conan build . --build=missing

Like the CMake build, this creates the executable coflasy-c.exe inside the build directory (Conan's intermediate build files are kept in build/Release).



## Exemplary Usage

The example.ini file contains all initial conditions and its structure should not be changed. Its entry `path` sets the location of the COFLASY main folder from which the input data (data/) is read and to which the results (output/) are written. A relative path is interpreted with respect to the directory the executable is run from. This means that the default `path = ..` works when running from the build directory
via

    ./coflasy-c.exe ../example.ini

For more explanation please see the file 'example.ini'.

## Contents

As of 14/10/2026, the GitHub repository of this code (https://github.com/mariofnavarro/COFLASY/tree/COFLASY-C) consists of the following source files and folders:

coflasy-c.cpp: contains the main body of the code.

source: contains auxiliary code and functions e.g. to read ini files, interpolations...

data: contains input data used by the code.

include: contains auxiliary libraries.

build: contains the .exe when the code is compiled.

CMakeLists: used to compile the code.

conanfile.py: used to build the code with Conan, which also fetches the GSL dependency.

example.ini: contains all initial conditions, which can be changed by the user.

output: stores the output of the code and a Mathematica notebook to generate simple plots.

publicdata: contains the data used in our paper [arXiv:2510.02438]


## Acknowledgement 

The code makes use of the ini file decoder [inifile-cpp](https://github.com/Rookfighter/inifile-cpp).

We thank Luca Wiggering (https://github.com/lucawiggering) for contributing the Conan installation method.

For any comment, please contact us at: Stefan.Sandner@lanl.gov, mario.fernandeznavarro@physik.uzh.ch

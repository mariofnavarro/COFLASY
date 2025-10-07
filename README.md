# COFLASY-C: COsmological evolution of FLAvour ASYmmetries - C++

The C++ code "COFLASY-C" has been developed by Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner in order to study the evolution of primordial lepton flavour asymmetries from the early Universe to BBN, following the approach of arXiv:2502.14960. If you use this code for scientific publications, please cite these papers: 

"Lepton Flavour Asymmetries: from the early Universe to BBN."  Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner,
[JHEP 06, 137] (https://link.springer.com/article/10.1007/JHEP06(2025)137), [arXiv:2502.14960] (https://arxiv.org/abs/2502.14960), [INSPIRE] (https://inspirehep.net/literature/2893306).

"A Limit on the Total Lepton Number in the Universe from BBN and the CMB." Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner, [arxiv:2510.02438]. (https://arxiv.org/abs/2510.02438).

## Requirements

| Package | Description |
| ------ | ----------- |
| [CMake](https://cmake.org)   | generate the Makefile |
| [GSL](https://www.gnu.org/software/gsl/) | numerical library for C++ |


## Installation

The Makefile can be generate within the main folder by executing the following:

    cd build
    cmake ..
    make

This should create the executable coflasy-c.exe inside the build directory.


## Exemplary Usage

    ./coflasy-c.exe ../example.ini

where the example.ini file contains all initial conditions and its structure should not be changed.
For more explanation please see the file 'example.ini'. 


## Contents

As of 03/10/2025, the GitHub repository of this code (https://github.com/mariofnavarro/COFLASY/tree/COFLASY-C) consists of the following source files and folders:

coflasy-c.cpp: contains the main body of the code.

source: contains auxiliary code and functions e.g. to read ini files, interpolations...

data: contains input data used by the code.

include: contains auxiliary libraries.

build: contains the .exe when the code is compiled.

CMakeLists: used to compile the code.

example.ini: contains all initial conditions, which can be changed by the user.

output: stores the output of the code and a Mathematica notebook to generate simple plots.

publicdata: contains the data used in our paper [arxiv:2510.XXXXX]


## Acknowledgement 

The code makes use of the ini file decoder [inifile-cpp](https://github.com/Rookfighter/inifile-cpp).


For any comment, please contact us at: Stefan.Sandner@lanl.gov, Mario.FernandezNavarro@glasgow.ac.uk

cmake_minimum_required(VERSION 3.12)

project(coflasy-c VERSION 1.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)


set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc-14" CACHE STRING "")
set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-14" CACHE STRING "")



find_package(Threads REQUIRED)
find_package(GSL REQUIRED)



include_directories(${GSL_INCLUDE_DIR})


add_executable( coflasy-c.exe ${CMAKE_CURRENT_SOURCE_DIR}/coflasy-c.cpp)
target_link_libraries(coflasy-c.exe GSL::gsl GSL::gslcblas)




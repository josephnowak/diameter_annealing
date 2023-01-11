cmake_minimum_required(VERSION 3.23)

set(CMAKE_CXX_STANDARD 17)

# set the project name
project(diameter_annealing)


# add the executable
add_library(diameter_annealing STATIC test_imports.cpp)


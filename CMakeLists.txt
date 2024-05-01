project(bhc3d LANGUAGES CXX)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")
set(CMAKE_CXX_FLAGS "-O3 -Wall -fopenmp -std=c++11 -fpermissive")

find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  message("${EIGEN3_INCLUDE_DIR}") 
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

find_package(LAPACK REQUIRED)

add_executable(ex_bhc3d main.cpp)
add_executable(ex_somig main_somig.cpp src/somig.cpp)
target_link_libraries(ex_somig
   LAPACK::LAPACK)

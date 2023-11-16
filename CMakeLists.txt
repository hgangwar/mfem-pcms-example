project(mfem-cmake-example C CXX)
cmake_minimum_required(VERSION 3.20)

find_package(MFEM REQUIRED)
find_package(pcms REQUIRED)
add_executable(fluxSolver fluxSolver.cpp)
target_link_libraries(fluxSolver PUBLIC mfem pcms::pcms)

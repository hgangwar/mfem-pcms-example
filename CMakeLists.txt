project(mfem-cmake-example C CXX)
cmake_minimum_required(VERSION 3.20)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(MFEM REQUIRED)
find_package(pcms REQUIRED)
find_package(spdlog REQUIRED)

add_executable(fluxSolver fluxSolver.cpp)
target_link_libraries(fluxSolver PUBLIC mfem pcms::pcms spdlog::spdlog)

add_executable(thermalSolver thermalSolver.cpp)
target_link_libraries(thermalSolver PUBLIC mfem pcms::pcms spdlog::spdlog)

add_executable(coupler mfem_pcms_coupling.cpp mfem_field_adapter.h)
target_link_libraries(coupler PUBLIC mfem pcms::pcms spdlog::spdlog)

add_executable(unit_test test_mfem_adapter.cpp)
target_link_libraries(unit_test PUBLIC mfem pcms::pcms spdlog::spdlog)

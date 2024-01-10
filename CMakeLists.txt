project(mfem-cmake-example C CXX)
cmake_minimum_required(VERSION 3.20)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(MFEM REQUIRED)
find_package(pcms REQUIRED)
#find_package(Catch2 3 REQUIRED)
find_program(BASH_EXECUTABLE NAMES bash REQUIRED)

add_executable(fluxSolver fluxSolver.cpp)
target_link_libraries(fluxSolver PUBLIC mfem pcms::pcms)

add_executable(thermalSolver thermalSolver.cpp)
target_link_libraries(thermalSolver PUBLIC mfem pcms::pcms)

add_executable(coupler mfem_pcms_coupling.cpp mfem_field_adapter.h)
target_link_libraries(coupler PUBLIC mfem pcms::pcms)


add_executable(mfem_adapter_unit_test test_mfem_adapter.cpp)
target_link_libraries(mfem_adapter_unit_test PUBLIC mfem pcms::pcms)

# add test_mfem_coupling as a test
add_executable(test_mfem_coupling test_mfem_coupling.cpp)
target_link_libraries(test_mfem_coupling PUBLIC mfem pcms::pcms)

# add this to enable testing
enable_testing()
# add the mfem_adapter_unit_test as a test
add_test(
  NAME mfem_adapter_unit_test 
  COMMAND $<TARGET_FILE:mfem_adapter_unit_test>
)

# add the test_mfem_coupling as a test
#add_test(
#  NAME test_mfem_coupling 
#  COMMAND $<TARGET_FILE:test_mfem_coupling>
#)

# add the test_mfem_coupling as a test using the bash script
add_test(
  NAME test_mfem_coupling 
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_mfem_coupling_bash_commands.sh $<TARGET_FILE:test_mfem_coupling>
)

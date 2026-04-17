cmake_minimum_required(VERSION 3.20)
project(pcms_mfem_coupling C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Omega_h REQUIRED)
#find_package(meshfields REQUIRED)
find_package(MFEM REQUIRED)
find_package(pcms REQUIRED)
#find_package(Gmsh REQUIRED)
find_program(BASH_EXECUTABLE NAMES bash REQUIRED)


add_executable(fluxSolver fluxSolver.cpp)
target_link_libraries(fluxSolver PUBLIC mfem pcms::pcms)

add_executable(thermalSolver thermalSolver.cpp)
target_link_libraries(thermalSolver PUBLIC mfem pcms::pcms) 

add_executable(coupler mfem_pcms_coupling.cpp mfem_field_adapter.h)
target_link_libraries(coupler PUBLIC mfem pcms::pcms)


add_executable(mfem_adapter_unit_test test_mfem_adapter.cpp)
target_link_libraries(mfem_adapter_unit_test PUBLIC mfem pcms::pcms)

add_executable(mfem_adapter_mask_test tests/test_mfem_mask.cpp)
target_link_libraries(mfem_adapter_mask_test PUBLIC mfem pcms::pcms)

add_executable(oh_adapter_mask_test tests/test_oh_mask_with_adapter.cpp)
target_link_libraries(oh_adapter_mask_test PUBLIC mfem pcms::pcms)

add_executable(describe describe.cpp)
target_link_libraries(describe PUBLIC mfem pcms::pcms)
# add test_mfem_coupling as a test
#add_executable(test_mfem_coupling test_mfem_coupling.cpp test_support.cpp)
#target_link_libraries(test_mfem_coupling PUBLIC mfem pcms::pcms  ${Gmsh_LIBRARIES})

# add test_mfem_coupling as a test
#add_executable(convg_test test_coupling_convg.cpp test_support.cpp)
#target_link_libraries(convg_test PUBLIC mfem pcms::pcms)

#add_executable(test_GDI test_GDI.cpp)
#target_link_libraries(test_GDI PUBLIC mfem pcms::pcms)


#add_executable(compare_mem_mesh compare_mem_mesh.cpp)
#target_link_libraries(compare_mem_mesh PUBLIC mfem pcms::pcms)

#add_executable(compare_serial_mesh compare_serial_mesh.cpp)
#target_link_libraries(compare_serial_mesh PUBLIC mfem pcms::pcms)

#add_executable(write_osh write_osh.cpp)
#target_link_libraries(write_osh PUBLIC mfem pcms::pcms)

add_executable(mfem_solve mfem_solve.cpp)
target_link_libraries(mfem_solve PUBLIC mfem)

add_executable(schwarz_coupling schwarz_coupling.cpp)
target_link_libraries(schwarz_coupling PUBLIC mfem)

add_executable(test_tag_transfer tests/test_tag_transfer.cpp)
target_link_libraries(test_tag_transfer PUBLIC mfem pcms::pcms)


add_executable(thermo_coupling thermo_coupling.cpp)
target_link_libraries(thermo_coupling PUBLIC mfem pcms::pcms)

add_executable(shift_mesh shift_mesh.cpp)
target_link_libraries(shift_mesh PUBLIC mfem pcms::pcms)

# add this to enable testing
enable_testing()

# add the mfem_adapter_unit_test as a test
#add_test(
#  NAME mfem_adapter_unit_test
#  COMMAND $<TARGET_FILE:mfem_adapter_unit_test>
#)

# add the test_mfem_coupling as a test
#add_test(
#  NAME test_mfem_coupling 
#  COMMAND $<TARGET_FILE:test_mfem_coupling>
#)

# add the test_mfem_coupling as a test using the bash script
#add_test(
#  NAME test_mfem_coupling 
#  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_mfem_coupling_bash_commands.sh $<TARGET_FILE:test_mfem_coupling>
#)
find_package(MPI REQUIRED)

#add_test(
#  NAME test_mfem_coupling
#  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
#          ${CMAKE_CURRENT_BINARY_DIR}/test_mfem_coupling
#)


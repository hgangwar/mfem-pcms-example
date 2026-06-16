cmake_minimum_required(VERSION 3.20)
project(pcms_mfem_coupling C CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Omega_h REQUIRED)
#find_package(meshfields REQUIRED)
find_package(MFEM REQUIRED)
find_package(pcms REQUIRED)
find_package(SUNDIALS REQUIRED)

find_package(MPI REQUIRED)

#find_package(Gmsh REQUIRED)
find_program(BASH_EXECUTABLE NAMES bash REQUIRED)

################### Coupling support #########################
add_library(coupling_support coupling_support.cpp)
target_include_directories(coupling_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(coupling_support PUBLIC pcms::pcms mfem)

add_library(mfem_adapter include/mfem_adapter_layout.cpp)
target_link_libraries(mfem_adapter PUBLIC pcms::pcms mfem)


################### MFEM SOLVERS #########################
#add_executable(fluxSolver fluxSolver.cpp)
#target_link_libraries(fluxSolver PUBLIC mfem pcms::pcms)

#add_executable(thermalSolver thermalSolver.cpp)
#target_link_libraries(thermalSolver PUBLIC mfem pcms::pcms)

add_executable(coupler mfem_pcms_coupling.cpp include/mfem_field_adapter.h)
target_link_libraries(coupler PUBLIC mfem pcms::pcms)

################### MFEM ADAPTER #########################
add_executable(mfem_adapter_unit_test test_mfem_adapter.cpp)
target_link_libraries(mfem_adapter_unit_test PUBLIC mfem pcms::pcms)

add_executable(mfem_adapter_mask_test tests/test_mfem_mask.cpp)
target_link_libraries(mfem_adapter_mask_test PUBLIC mfem pcms::pcms)

################### OH ADAPTER TEST #########################
add_executable(oh_adapter_mask_test tests/test_oh_mask_with_adapter.cpp)
target_link_libraries(oh_adapter_mask_test PUBLIC mfem pcms::pcms)

add_executable(test_permutation tests/test_permutation.cpp)
target_link_libraries(test_permutation PRIVATE coupling_support pcms::pcms mfem)

################### Mesh Operations #########################
add_executable(describe utilities/describe.cpp)
target_link_libraries(describe PUBLIC mfem pcms::pcms)

add_executable(shift_mesh utilities/shift_mesh.cpp)
target_link_libraries(shift_mesh PUBLIC mfem pcms::pcms)

add_executable(test_tag_transfer tests/test_tag_transfer.cpp)
target_link_libraries(test_tag_transfer PUBLIC mfem pcms::pcms)

################### Schwartz Iterations #########################

add_executable(schwarz_coupling tests/schwarz_coupling.cpp)
target_link_libraries(schwarz_coupling PUBLIC mfem)


add_executable(thermo_coupling thermo_coupling.cpp)
target_link_libraries(thermo_coupling PUBLIC mfem pcms::pcms coupling_support)

# Full overlap, no mask coupling convergence test
add_executable(coup_convg tests/test_coupling_convg.cpp)
target_link_libraries(coup_convg PUBLIC mfem pcms::pcms)


################### SUNDIALS Schwartz Iterations #########################
add_library(schwarz_sundial Schwarz_Sundial_Coupling.cpp)

target_include_directories(schwarz_sundial PUBLIC
        ${SUNDIALS_ROOT}/include
)

target_link_libraries(schwarz_sundial PUBLIC
        mfem
        ${SUNDIALS_ROOT}/lib64/libsundials_core.so
        ${SUNDIALS_ROOT}/lib64/libsundials_nvecserial.so
        ${SUNDIALS_ROOT}/lib64/libsundials_arkode.so
)
add_executable(stepper Sunstepper.cpp)

target_link_libraries(stepper PUBLIC schwarz_sundial)

add_executable(pcms_sun_driver
        test_pcms_sundial.cpp
)

target_link_libraries(pcms_sun_driver PUBLIC
        mfem
        pcms::pcms
        schwarz_sundial
)





export NVCC_WRAPPER_DEFAULT_COMPILER=`which mpicxx`
cmake -S . -B build \
	-DCMAKE_CXX_COMPILER=`which mpicxx` \
	-DCMAKE_C_COMPILER=`which mpicc` \
	-DCMAKE_INSTALL_PREFIX=$HOME/src/mfem-pcms-example/install \
	-DMFEM_ROOT=$HOME/src/MFEM/mfem/install \
	-Dpcms_ROOT=$HOME/src/PCMS/build/ADA89/pcms/install \
	-DOmega_h_ROOT=$HOME/src/PCMS/deps/build/ADA89/omega_h/install/ \
	-Dperfstubs_DIR=$HOME/src/PCMS/deps/build/perfstubs/install/lib/cmake/ \
	-DCMAKE_PREFIX_PATH=/users/gangwh/src/gmsh/build/install/ \
	-Dmeshfields_ROOT=/users/gangwh/src/PCMS/deps/build/ADA89/meshFields/install/ \
	-DKokkos_ROOT=/users/gangwh/src/PCMS/deps/build/ADA89/kokkos/install \
	-DCMAKE_BUILD_TYPE=Debug

  	#-DGmsh_LIBRARIRES=/users/gangwh/src/gmsh/build/install/lib64/ \
  	#-DGmsh_INCLUDE_DIRS=/users/gangwh/src/gmsh/build/install/include \


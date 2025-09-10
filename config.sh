export NVCC_WRAPPER_DEFAULT_COMPILER=`which mpicxx`
cmake -S . -B build \
	-DCMAKE_CXX_COMPILER=`which mpicxx` \
	-DCMAKE_C_COMPILER=`which mpicc` \
	-DCMAKE_INSTALL_PREFIX=$HOME/src/mfem-pcms-example/install \
	-DMFEM_ROOT=$HOME/src/MFEM/mfem/install \
	-Dpcms_ROOT=$HOME/src/PCMS/build/ADA89/pcms/install \
	-Dspdlog_ROOT=$HOME/src/PCMS/deps/build/spdlog/install \
	-Dperfstubs_DIR=$HOME/src/PCMS/deps/build/perfstubs/install/lib/cmake/ \
  	-DCMAKE_BUILD_TYPE=RelWithDebInfo \

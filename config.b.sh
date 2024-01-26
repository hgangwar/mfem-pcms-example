#export NVCC_WRAPPER_DEFAULT_COMPILER=`which mpicxx`
cmake -S . -B build \
	-DCMAKE_CXX_COMPILER=`which mpicxx` \
	-DCMAKE_C_COMPILER=`which mpicc` \
	-DMFEM_ROOT=/lore/hasanm4/mFEM/install/mfem-mpi/ \
	-Dpcms_ROOT=/lore/hasanm4/pcms-coupler/build/OMP/pcms/install \
	-Dperfstubs_DIR=/lore/mersoj/laces-software/build/perfstubs/install/lib64/cmake/ \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DGmsh_ROOT=/lore/hasanm4/Gmsh/ \
#  -DCatch2_ROOT=/lore/hasanm4/Catch2/ \
#	-Dperfstubs_ROOT=/lore/mersoj/laces-software/build/perfstubs/install

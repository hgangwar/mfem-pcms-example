project=/users/gangwh/src/mfem-pcms-example
cd $project/build
rm -rf *.bp
rm -rf *.sol*
echo "removed old files from current directory"

echo "Running the test."

mpirun -np 1 ./mfem_only_test 1  $project/mesh/cube.msh CG Jacobi & \
mpirun -np 1 ./mfem_only_test -1 $project/mesh/cube.osh  & \
mpirun -np 1 ./mfem_only_test 0  $project/mesh/cube.msh CG Jacobi

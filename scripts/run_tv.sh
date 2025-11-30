project=/users/gangwh/src/mfem-pcms-example
cd $project/build
rm -rf *.bp
echo "removed old bp files"
echo "Running the test."

totalview --args \
mpirun -np 1 ./convg_test 1  $project/mesh/cube.msh CG Jacobi & \
mpirun -np 1 ./convg_test -1 $project/mesh/cube.osh  & \
mpirun -np 1 ./convg_test 0  $project/mesh/cube.msh CG Jacobi


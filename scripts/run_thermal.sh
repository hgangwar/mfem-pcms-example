project=/users/gangwh/src/mfem-pcms-example
cd $project/build
rm -rf *.bp
echo "removed old bp files"

echo "Running the test."

mpirun -np 1 ./tt_coupling 1  $project/mesh/cube.msh CG Jacobi & \
mpirun -np 1 ./tt_coupling -1 $project/mesh/cube.osh  & \
mpirun -np 1 ./tt_coupling 0  $project/mesh/cube.msh CG Jacobi

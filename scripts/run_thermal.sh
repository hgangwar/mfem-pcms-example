project=/users/gangwh/src/mfem-pcms-example
cd $project/build
rm -rf *.bp
echo "removed old bp files"
echo "Running the test."

mpirun -np 1 ./thermo_coupling 0  $project/mesh/box_tri.msh CG Jacobi & \
mpirun -np 1 ./thermo_coupling -1 $project/mesh/box_tri.osh  & \
mpirun -np 1 ./thermo_coupling 1  $project/mesh/box_tri.msh CG Jacobi

project=/users/gangwh/src/mfem-pcms-example
cd $project/build
echo $project
rm -rf *.bp
echo "removed old bp files"
echo "Running the test."


mpirun -np 1 ./thermo_coupling 0  $project/mesh/box_A_marked.mesh CG Jacobi & \
mpirun -np 1 ./thermo_coupling -1 $project/mesh/box_A_marked.osh $project/mesh/box_B_marked.osh & \
mpirun -np 1 ./thermo_coupling 1  $project/mesh/box_B_marked.mesh CG Jacobi

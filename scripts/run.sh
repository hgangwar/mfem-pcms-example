project=/users/gangwh/src/mfem-pcms-example
cd $project/build
rm -rf *.bp
echo "removed old bp files"

echo "Running the test."

mpirun -np 1 ./coup_convg 1  $project/mesh/box_A.msh CG Jacobi & \
mpirun -np 1 ./coup_convg -1 $project/mesh/box_A.osh  & \
mpirun -np 1 ./coup_convg 0  $project/mesh/box_A.msh CG Jacobi

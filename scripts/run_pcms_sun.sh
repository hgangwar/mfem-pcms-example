project=/users/gangwh/src/mfem-pcms-example
cd $project/build
echo $project
rm -rf *.bp
echo "removed old bp files"
echo "Running the test."


mpirun -np 1 ./pcms_sun_driver 1  $project/mesh/box_A.msh CG Jacobi & \
mpirun -np 1 ./pcms_sun_driver -1 $project/mesh/box_A.osh  & \
mpirun -np 1 ./pcms_sun_driver 0  $project/mesh/box_A.msh CG Jacobi

echo "removed old bp files"
rm -rf *.bp
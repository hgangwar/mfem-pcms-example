project=/users/gangwh/src/mfem-pcms-example
cd $project/build || exit
echo $project

rm -rf *.bp
echo "removed old bp files"
echo "Running the test."


totalview -args mpirun -np 1 ./pcms_sun_driver 1  $project/mesh/box_A.msh CG Jacobi & \
totalview -args mpirun -np 1 ./pcms_sun_driver -1 $project/mesh/box_A.osh  & \
totalview -args mpirun -np 1 ./pcms_sun_driver 0  $project/mesh/box_A.msh CG Jacobi

echo "removed old bp files"
rm -rf *.bp
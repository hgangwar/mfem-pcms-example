rm -rf *.bp
echo "removed old bp files"
echo "Running the test."
totalview --args mpirun \
-np 1 ./build/convg_test 0  mesh/cube.msh  : \
-np 1 ./build/convg_test -1 mesh/cube.osh  : \
-np 1 ./build/convg_test 1  mesh/cube.msh 


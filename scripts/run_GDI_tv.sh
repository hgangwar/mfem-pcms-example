cd $HOME/src/mfem-pcms-example/build
rm -rf *.bp
echo "removed old bp files"
echo "Running the test."

totalview --args mpirun \
-np 1 ./test_GDI 1 : \
-np 1 ./test_GDI -1 : \
-np 1 ./test_GDI 0  
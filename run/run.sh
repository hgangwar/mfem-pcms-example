# Go to ..
cd ..

# build the project
# rm -r build
source ./config.sh

cd build
make -j 8

# run the project
cd ../run

#totalview mpirun -a -np 4 ../build/fluxSolver
mpirun -np 4 ../build/fluxSolver

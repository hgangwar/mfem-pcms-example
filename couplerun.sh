mkdir Run
cd Run

mpirun -np 1 ../build/fluxSolver &
mpirun -np 1 ../build/thermalSolver &
mpirun -np 1 ../build/coupler

cd ..

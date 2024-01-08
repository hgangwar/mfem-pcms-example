#! /bin/bash
# delete if Run directory already exists
#if [ -d "Run" ]; then
#  rm -rf Run
#fi

# create run directory if it doesn't exist
if [ ! -d "Run" ]; then
  mkdir Run
fi
cd Run

declare -a PIDS=()

kill_stuff() {
  for PID in "${PIDS[@]}"
  do
    echo "Killing $PID"
    kill -9 $PID
  done
  cd ..
}


run_stuff() {
  mpirun -np 1 $1 &
  PIDS+=($!)
}

trap "kill_stuff" SIGINT

run_stuff ../build/fluxSolver &> fluxSolver.log &
run_stuff ../build/thermalSolver &> thermalSolver.log &
#run_stuff ../build/coupler &> coupler.log

#run_stuff ../build/thermalSolver &
#run_stuff ../build/fluxSolver &
#totalview -args mpirun  -np 1 ../build/coupler 
run_stuff ../build/coupler 


wait



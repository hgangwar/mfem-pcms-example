#! /bin/bash
# delete if Run directory already exists
if [-d Run]; then rm -r Run; fi
mkdir Run
cd Run

declare -a PIDS=()

kill_stuff() {
  for PID in "${PIDS[@]}"
  do
    echo "Killing $PID"
    kill -9 $PID
  done
}


run_stuff() {
  mpirun -np 1 $1 &
  PIDS+=($!)
}

trap "kill_stuff" SIGINT

run_stuff ../build/fluxSolver
run_stuff ../build/thermalSolver
run_stuff ../build/coupler

wait


##mpirun -np 1 ../build/fluxSolver &
#PIDS += $!
#mpirun -np 1 ../build/thermalSolver &
#PIDS += $!
#mpirun -np 1 ../build/coupler&
#PIDS += $!
#
#for PID in "$PIDS[@]" do
#  kill -9 $PID
#done
#
#wait
#
## return back to the case home
#cd ..

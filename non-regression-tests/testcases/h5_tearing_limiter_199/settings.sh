#!/bin/bash

description="Test case for time evolution of model 199: Tearing mode in circular plasma at n_tor=3."

# --- Model used
jorekmodel="199"

# --- Files required to run the code (executables copied automatically)
requiredfiles="$codedir/jorek_model${jorekmodel}_1 $codedir/jorek_model${jorekmodel}_3"

# --- How many MPI tasks and OpenMP threads are required?
mpitasks=2
ompthreads=4

# --- Compare which kind of data? (separate keywords with blanks)
comparedata="???"

# --- Tolerance for the comparison. A comparison fails if both thresholds are exceeded.
relativeaccuracy="???"
absoluteaccuracy=""

function compile_jorek () {
    ./util/config.sh model=$jorekmodel "n_tor=1 n_plane=1 n_period=1"
    make clean && make -j 3 jorek_model${jorekmodel} &&\
    cp jorek_model${jorekmodel} jorek_model${jorekmodel}_1
    returncode=$?
    if [ $returncode -eq 0 ]; then
	./util/config.sh model=$jorekmodel "n_tor=3 n_plane=4 n_period=1"
	make clean &&  make -j 3 jorek_model${jorekmodel} &&\
      cp jorek_model${jorekmodel} jorek_model${jorekmodel}_3  
	returncode=$?
    fi
    if [ $returncode -eq 0 ]; then
	./util/config.sh model=$jorekmodel "n_tor=3 n_plane=4 n_period=1 "
	make clean; 
	returncode=$?
    fi
    return $returncode
}

function run_jorek_first () {
  sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;"  $codedir/namelist/model199/intear > input
  eval $PRERUN
  export OMP_NUM_THREADS=$ompthreads
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=0 tstep_n=1 n_flux=35 n_tht=14
  $MPIRUN 1 ./jorek_model${jorekmodel}_1 < input > logfile &&\
    cp jorek_restart.rst jorek_equil.rst || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=89 tstep_n=1000 n_flux=35 n_tht=14
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input > logfile &&\
    cp jorek00089.rst jorek00089_export.rst || exit 1
  
}

function run_jorek_second () {
  sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;"  $codedir/namelist/model199/intear > input
  eval $PRERUN
  export OMP_NUM_THREADS=$ompthreads
  cp ${testcasedir}/jorek00089_export.rst jorek_restart.rst || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1000 n_flux=35 n_tht=14
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input > logfile &&\
    cp jorek00090.rst jorek00090_export.rst || exit 1
  
}

function run_jorek () {
#    run_jorek_first;
    run_jorek_second;
}

function compare_jorek_res () {
  returncode=0
  return $returncode
}
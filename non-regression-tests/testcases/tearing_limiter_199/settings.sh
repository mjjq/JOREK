#!/bin/bash

description="Test case for time evolution of model 199: Tearing mode in circular plasma at n_tor=3."

# --- Model used
jorekmodel="199"

# --- Files required to run the code (executables copied automatically)
requiredfiles="$codedir/jorek_model${jorekmodel}_1 $codedir/jorek_model${jorekmodel}_3 $codedir/rst_bin2hdf5"

# --- How many MPI tasks and OpenMP threads are required?
mpitasks=2
ompthreads=4

function compile_jorek () {
    returncode=0
    ./util/config.sh model=$jorekmodel "n_tor=1 n_plane=1 n_period=1"
    make clean && make -j 3 jorek_model${jorekmodel} &&\
    cp jorek_model${jorekmodel} jorek_model${jorekmodel}_1
    returncode=$?
    if [ $returncode -eq 0 ]; then
	./util/config.sh model=$jorekmodel "n_tor=3 n_plane=4 n_period=1"
	make clean &&  make -j 3 jorek_model${jorekmodel} &&\
          make  rst_bin2hdf5 &&\
          cp jorek_model${jorekmodel} jorek_model${jorekmodel}_3  
	returncode=$?
    fi
    return $returncode
}


function restart_run () {
  sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;"  $codedir/namelist/model199/intear > input
  if [ -n "$PRERUN" ]; then
    eval $PRERUN
  fi
  export OMP_NUM_THREADS=$ompthreads
  cp ${testcasedir}/jorek00089_export.rst jorek_restart.rst || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1000 n_flux=35 n_tht=14
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile   
}

function initial_run () {
  sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;"  $codedir/namelist/model199/intear > input
  if [ -n "$PRERUN" ]; then
    eval $PRERUN
  fi
  export OMP_NUM_THREADS=$ompthreads
  # Equilibrium computation
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=0 tstep_n=1 n_flux=35 n_tht=14
  ($MPIRUN 1 ./jorek_model${jorekmodel}_1 < input > logfile &&\
    cp jorek_restart.rst jorek_equil.rst) || exit 1
  # Time evolution
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=89 tstep_n=1000 n_flux=35 n_tht=14
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile || exit 1
  cp jorek00089.rst ${testcasedir}/jorek00089_export.rst || exit 1
  # Restart run (1 time step)
  restart_run
  # Export the final restart file as an HDF5 file
  tstep=90;
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  cp jorek00090.h5 ${testcasedir}/jorek00090_export.h5 || exit 1
}

function compare_results () {
  # convert binary restart file into hdf5 file
  tstep=90;
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  # compare with reference file and return the result
  h5diff -d 1e-15 jorek00090.h5 ${testcasedir}/jorek00090_export.h5 values; 
  returncode=$?
  return $returncode
}

function pack_restart_files () {
  cd ${testcasedir} || exit 1  
  testname=$(basename $testcasedir)
  tar cvzf ${testname}.tgz jorek00090_export.h5 jorek00089_export.rst || exit 1
}
#!/bin/bash

description="Time evolution of model 199: Tearing mode in circular plasma at n_tor=3."

# --- Model used
jorekmodel="199"

# --- Files required to run the code (executables copied automatically)
requiredfiles="$testcasedir/jorek_model${jorekmodel}_1 $testcasedir/jorek_model${jorekmodel}_3 $testcasedir/rst_bin2hdf5 $testcasedir/rst_hdf52bin $testcasedir/input"

# --- How many MPI tasks and OpenMP threads are required?
mpitasks=2
ompthreads=16

function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1
  make cleanall                                             || exit 1
  make $compilopt                                           || exit 1
  cp jorek_model${jorekmodel} jorek_model${jorekmodel}_1    || exit 1
  
  ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=1
  make clean                                                || exit 1
  make $compilopt                                           || exit 1
  make $compilopt rst_bin2hdf5 rst_hdf52bin                 || exit 1
  cp jorek_model${jorekmodel} jorek_model${jorekmodel}_3    || exit 1
  
  cp jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3 rst_hdf52bin rst_bin2hdf5 $testcasedir/ || exit 1
}

function restart_run () {
  if [ -n "$PRERUN" ]; then
    eval $PRERUN
  fi
  export OMP_NUM_THREADS=$ompthreads
  
  ### REMARK MHOELZL: THIS PART SHOULD GO INTO A BASH FUNCTION IN run_test.sh
  ### SUCH THAT WE ONLY CALL THE FOLLOWING
  ### use_restartfile 00089
  cp ${testcasedir}/jorek00089_export.h5 jorek00089.h5 || exit 1
  tstep=89
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_hdf52bin
  rm file.out
  cp jorek00089.rst jorek_restart.rst
  ### END SECTION
  
  ### REMARK MHOELZL: PARAMETERS LIKE n_flux, n_tht etc should simply be set
  ### in the input file, not here
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1000
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile   
}

function initial_run () {
  if [ -n "$PRERUN" ]; then
    eval $PRERUN
  fi
  export OMP_NUM_THREADS=$ompthreads

  # Equilibrium computation
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=0
  ($MPIRUN 1 ./jorek_model${jorekmodel}_1 < input > logfile &&\
    cp jorek_restart.rst jorek_equil.rst) || exit 1

  # Time evolution
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=89 tstep_n=1000
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile || exit 1

  ### REMARK MHOELZL: THE FOLLOWING SHOULD ALSO BE ENCAPSULATED IN A BASH FUNCTION SUCH
  ### THAT WE ONLY CALL SOMETHING LIKE:
  ### export_restart 00089
  # Export final restart file as HDF5 file
  tstep=89
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  cp jorek000${tstep}.h5 ${testcasedir}/jorek000${tstep}_export.h5 || exit 1

  # Restart run (1 time step)
  restart_run

  ### REMARK MHOELZL: AGAIN ENCAPSULATE IN BASH FUNCTION
  # Export  final restart file as HDF5 file
  tstep=90
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  cp jorek000${tstep}.h5 ${testcasedir}/jorek000${tstep}_export.h5 || exit 1
}

function compare_results () {
  # convert binary restart file into hdf5 file
  tstep=90;
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  # compare with reference file and return the result
  h5diff -d 1e-14 jorek00090.h5 ${testcasedir}/jorek00090_export.h5 values; 
  returncode=$?
  return $returncode
}

function pack_restart_files () {
  cd ${testcasedir} || exit 1  
  testname=$(basename $testcasedir)
  tar cvzf ${testname}.tgz jorek00090_export.h5 jorek00089_export.h5 || exit 1
}

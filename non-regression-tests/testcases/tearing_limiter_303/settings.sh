#!/bin/bash

echo " IN testcases/tearing_limiter_303/settings.sh "
echo " ============================================ "

description="Test case for time evolution of model 303: Tearing mode in circular plasma at n_tor=3."

# --- Model used
jorekmodel="303"
COMMONOPT="n_flux=26 n_tht=28 n_open=8 n_leg=8 n_private=6"

# --- Files required to run the code (executables copied automatically)
requiredfiles="$codedir/jorek_model${jorekmodel}_1 $codedir/jorek_model${jorekmodel}_3 $codedir/rst_bin2hdf5 $testcasedir/input"

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
          make rst_bin2hdf5 &&\
          cp jorek_model${jorekmodel} jorek_model${jorekmodel}_3  
	returncode=$?
    fi
    return $returncode
}


function restart_run () {
  sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;"  $testcasedir/input > input
  if [ -n "$PRERUN" ]; then
    eval $PRERUN
  fi
  export OMP_NUM_THREADS=$ompthreads
  cp ${testcasedir}/jorek00940_export.rst jorek_restart.rst || exit 1
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=10 $COMMONOPT nout=1
  ($MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_end &&\
   cp jorek00941.rst jorek00941_export.rst|| exit 1)
}

function initial_run () {
  sed "s/nstep.*=/nstep_n =/;s/tstep.*=/tstep_n =/;"  $testcasedir/input > input
  if [ -n "$PRERUN" ]; then
    eval $PRERUN
  fi
  export OMP_NUM_THREADS=$ompthreads

  echo " mpitasks = " $mpitasks
  # Equilibrium computation
# => jorek_equil.rst
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=0 tstep_n=1 $COMMONOPT
  cat input
  echo " we are HERE !!" $PWD
  echo " mpirun " $MPIRUN
  echo " mpitasks " $mpitasks
  ($MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input > logfile_equil &&\
    cp jorek_restart.rst jorek_equil.rst) || exit 1

  # Time evolution
  # => jorek00120.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 10, 10, 10, 10, 40, 40' 'tstep_n= 1e-3, 1e-2, 1e-1, 5e-1, 1, 2' $COMMONOPT nout=10
  ($MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input >> logfile_1 &&\
    cp jorek_restart.rst jorek_rst1.rst) || exit 1
  # => jorek00240.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 20, 100' 'tstep_n= 2, 5' $COMMONOPT  nout=10
  ($MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_2 &&\
    cp jorek_restart.rst jorek_rst2.rst) || exit 1
  # => jorek00340.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_3 || exit 1
  # => jorek00440.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_4 || exit 1
  # => jorek00540.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_4 || exit 1
  # => jorek00640.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_5 || exit 1
  # => jorek00740.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_6 || exit 1
  # => jorek00840.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_7 || exit 1
  # => jorek00940.rst
  ${codedir}/util/setinput.sh input restart=.t. 'nstep_n= 100' 'tstep_n= 10' $COMMONOPT  nout=10
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_8 || exit 1
  # => jorek00941.rst
  restart_run

  # Export the final restart file as an HDF5 file
  tstep=941
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  cp jorek00941.h5 ${testcasedir}/jorek00941_export.h5 || exit 1
}

function compare_results () {
  # convert binary restart file into hdf5 file
  printf "%5.5d \n" $tstep > ./file.out
  ./rst_bin2hdf5
  rm file.out
  # compare with reference file and return the result
  h5diff -d 1e-15 jorek00941.h5 ${testcasedir}/jorek00941_export.h5 values; 
  returncode=$?
  return $returncode
}

function pack_restart_files () {
  cd ${testcasedir} || exit 1  
  testname=$(basename $testcasedir)
  tar cvzf ${testname}.tgz jorek00941_export.h5 jorek00940_export.rst || exit 1
}
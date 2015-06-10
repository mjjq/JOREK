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
    make $compilopt  jorek_model${jorekmodel}                 || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1    || exit 1
    
    ./util/config.sh model=$jorekmodel n_tor=3 n_plane=4 n_period=1
    make cleanall                                             || exit 1
    make $compilopt jorek_model${jorekmodel} \
	rst_bin2hdf5 rst_hdf52bin                             || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3    || exit 1
    
    cp jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3 rst_hdf52bin rst_bin2hdf5 $testcasedir/ || exit 1
}

function restart_run () {
    if [ -n "$PRERUN" ]; then
	eval $PRERUN
    fi
    export OMP_NUM_THREADS=$ompthreads
    
    # Import restart file
    cp ${testcasedir}/jorek00089_export.h5 jorek_restart.h5 || exit 1
    ./rst_hdf52bin # generate jorek_restart.rst
 
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
    
    # Export final restart file as HDF5 file
    ./rst_bin2hdf5 # take the jorek_restart.rst as input
    cp jorek_restart.h5 ${testcasedir}/jorek00089_export.h5 || exit 1

    # Restart run (1 time step)
    restart_run

    # Export final restart file as HDF5 file
    ./rst_bin2hdf5 # take the jorek_restart.rst as input
    cp jorek_restart.h5 ${testcasedir}/jorek00090_export.h5 || exit 1
}

function compare_results () {
    # convert last binary restart file into hdf5 file    
    ./rst_bin2hdf5  # jorek00090.rst taken as input

    # compare with reference file and return the result
    h5diff -d 1e-14 jorek_restart.h5 ${testcasedir}/jorek00090_export.h5 values; 
    returncode=$?
    return $returncode
}

function pack_restart_files () {
    cd ${testcasedir} || exit 1  
    testname=$(basename $testcasedir)
    tar cvzf ${testname}.tgz jorek00090_export.h5 jorek00089_export.h5 || exit 1
}


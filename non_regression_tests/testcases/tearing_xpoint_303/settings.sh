#!/bin/bash

description="Time evolution of model 303: Tearing mode in circular plasma at n_tor=3."
jorekmodel="303"
mpitasks=2
ompthreads=16
requiredfiles="jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3 rst_bin2hdf5 rst_hdf52bin input"


# --- Compile the code for the test case
function compile_jorek () {
    ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1
    make cleanall                                             || exit 1
    make $compilopt jorek_model${jorekmodel}                  || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1    || exit 1

    ./util/config.sh model=$jorekmodel n_tor=3 n_plane=8 n_period=2
    make cleanall                                             || exit 1
    make $compilopt jorek_model${jorekmodel} \
	rst_bin2hdf5 rst_hdf52bin                             || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_3    || exit 1
    
    cp jorek_model${jorekmodel}_1 jorek_model${jorekmodel}_3 rst_hdf52bin rst_bin2hdf5 $testcasedir/ || exit 1
}

# --- Re-run the whole case from scratch into the non-linear phase
function initial_run () {
    # Equilibrium computation
    # => jorek_equil.rst
    ${codedir}/util/setinput.sh input restart=.f. nstep_n=0 tstep_n=1
    $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input > logfile_equil || exit 1
    cp jorek_restart.rst jorek_equil.rst || exit 1

    # Time evolution
    # => jorek00120.rst
    ${codedir}/util/setinput.sh input restart=.t. nstep_n=10, 10, 10, 10, 40, 40 tstep_n=1e-3, 1e-2, 1e-1, 5e-1, 1, 2 nout=10
    $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_1 || exit 1
    cp jorek_restart.rst jorek_rst1.rst || exit 1
    # => jorek00240.rst
    ${codedir}/util/setinput.sh input restart=.t. nstep_n=20, 100 tstep_n=2, 5 nout=10
    $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_2 || exit 1
    cp jorek_restart.rst jorek_rst2.rst || exit 1
    # => jorek00940.rst
    ${codedir}/util/setinput.sh input restart=.t. nstep_n=700 tstep_n=10 nout=10
    $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile_3 || exit 1
}

# --- Carry out the test case, i.e., run a single time step in the non-linear phase
function restart_run () {
    ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=10  || exit 1
    $MPIRUN $mpitasks ./jorek_model${jorekmodel}_3 < input >> logfile   || exit 1
}

# --- Compare the results of the test case to the reference solution
function compare_results () {
    ./rst_bin2hdf5                                                      || exit 1 
    h5diff -p 1e-6 jorek_restart.h5 ${testcasedir}/end.h5 values        || exit 1 
}


# --- General settings
jorekmodel="199"
description="VDE test case for an ITER plasma with simplified wall geometry (JOREK-STARWALL, model199)."
mpitasks=1
binaries="jorek_model${jorekmodel}_1 rst_bin2hdf5 rst_hdf52bin"
requiredfiles="$binaries input starwall-response.dat coil_field.dat"
extra_remote_files="starwall-response.dat coil_field.dat"


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1                    || exit 1
  make cleanall                                                                      || exit 1
  make $compilopt $debugoptions jorek_model${jorekmodel} rst_bin2hdf5 rst_hdf52bin   || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1
}


# --- Re-run the whole case from scratch
function initial_run () {
  ${codedir}/util/setinput.sh input restart=.f. nstep_n=5,5,10,10,90,10,170 tstep_n=1.d0,1.d1,1.d2,3.d2,1.d3,3.d2,1.d2 || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile               || exit 1
}


# --- Carry out the test case.
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.t. nstep_n=1 tstep_n=1.d2 || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee -a logfile            || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  ./rst_bin2hdf5 < ./input                                                           || exit 1
  h5diff -d 1e-5 jorek_restart.h5 ${testcasedir}/end.h5 values                      || exit 1
}

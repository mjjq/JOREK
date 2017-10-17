# --- General settings
jorekmodel="199"
description="Free boundary equilibrium for ASDEX Upgrade plasma (JOREK-STARWALL)."
mpitasks=2
binaries="jorek_model${jorekmodel}_1 rst_bin2hdf5 rst_hdf52bin"
binaries_initial=""
requiredfiles="input starwall-response.dat coil_field.dat ffp.dat t.dat rho.dat"
extra_remote_files="starwall-response.dat coil_field.dat"


# --- Compile the code for the test case
function compile_jorek () {
  ./util/config.sh model=$jorekmodel n_tor=1 n_plane=1 n_period=1                    || exit 1
  make cleanall                                                                      || exit 1
  make $compilopt jorek_model${jorekmodel} rst_bin2hdf5 rst_hdf52bin                 || exit 1
  mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                             || exit 1
}


# --- Re-run the whole case from scratch into the non-linear phase
function initial_run () {
  # --- Dummy run not actually needed...
  ${codedir}/util/setinput.sh input restart=.f. freeboundary_equil=.f.               || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile               || exit 1
}


# --- Carry out the test case (not actually a restart for this testcase...).
function restart_run () {
  ${codedir}/util/setinput.sh input restart=.f. freeboundary_equil=.t.               || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee -a logfile            || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  ./rst_bin2hdf5 < ./input                                                           || exit 1
  h5diff -d 1e-9 jorek_restart.h5 ${testcasedir}/end.h5 values                       || exit 1
}

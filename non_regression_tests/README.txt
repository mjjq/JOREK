
1) Load the needed module and put in your PATH the compiler you will use.
Copy/paste a Makefile.inc that works for the machine you are using in jorek directory.

2) Download the reference test cases (restart files) from web site
$ non-regression-tests/get-all-data.sh

3) Compile the jorek executables for the test "tearing_limiter_199" using the command:
$ non-regression-tests/run_test.sh -p tearing_limiter_199

4) Run the simulation using a batch script that calls:
$ non-regression-tests/run_test.sh -k -n tearing_limiter_199

   As an example, you can look at job_scripts/poincare/tearing199.job
   to design a batch script that does it on a parallel machine.

Optional) To send a new set of restart files to the web site that encompass
   the database :
$  cd non-regression-tests/testcase; ./send-testcase-data.sh tearing_limiter_199

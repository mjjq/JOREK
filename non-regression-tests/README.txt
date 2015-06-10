
Load the needed module and put in your PATH the compiler you will use.
Copy/paste a Makefile.inc that works for the machine you are using in jorek directory.

Compile the jorek executables using the command:
$ non-regression-tests/run_test.sh -r tearing_limiter_199

Run the simulation using a batch script that calls:
$ non-regression-tests/run_test.sh -n tearing_limiter_199

For example, you can look at examples/poincare/tearing199.cmd

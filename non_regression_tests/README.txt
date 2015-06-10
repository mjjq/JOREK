Non regression and regression tools for Jorek
---------------------------------------------

1) INTRODUCTION
===============

Tests should be devised to detect bugs that might lead to crashes as
early as possible so as to ease their correction. Software
verification and validation is a continual effort to ensure that the
simulation code meets all design requirements and that computed
results provide a high level of accuracy.

For Jorek code, a strategy has been setup to launch several automated
tests when required, typically at each commit on the software
repository. In addition to the automated tests, the users and
developpers have also the possibility to launch the same tests by
themselves to check their own version of the code.

2) IN A NUTSHELL
================

Here is a quick overview of how to perform the tests as a user.
As an example, the setup on the Helios machine is given.

1) Get the code from the repository with git 
   on a machine where outward ssh connection is available

      $ git clone ssh://git.iter.org/stab/jorek.git jorek_git
      $ git checkout -b nrt feature/nrt

2) Download the reference test cases (restart files) from the
   jorek.eu web site 

      $ jorek_git/non_regression_tests/get_all_data.sh
  
3) Transfer the directory to Helios 
     
      $ rsync -av --progress jorek_git USER@helios.iferc-csc.org:
      $ ssh USER@helios.iferc-csc.org
      $ cd jorek_git

4) Load the needed module and put in your PATH the compiler you will use.

      $ source non_regression_tests/job_scripts/helios/env.sh

   Copy/paste a Makefile.inc that works for the machine you are using in the
   jorek directory.
 
      $ cp Make.inc/Makefile.helios_h5_linux_intel.inc Makefile.inc

5) Compile the jorek executables for the test "tearing_xpoint_303" 
   using the command:

      $ non_regression_tests/run_test.sh -p tearing_xpoint_303

6) Run the simulation using a batch script that calls:

      $ non_regression_tests/run_test.sh -k -n tearing_limiter_199

   On Helios you can for example submit the following job

      $ cd non_regression_tests/job_scripts/
      $ sbatch tearing_xpoint_303.job
   
   If the test is successful the end of the output file should be

      $ tail -n 1 tearing_*out
       Test 'tearing_xpoint_303' passed.

3) DIRECTORIES AND FILES
========================

Here is a short description of the centents of the 
non_regression_test directory: 

 get_all_data.sh : download all reference restart files from jorek.eu
 README.txt      : this file
 run_test.sh     : wrapper script that can initiate compilation or
                   launch mpirun command on a testcase given as parameter
 compile_all.sh  : compile all testcases one after the other 
                   (run_test.sh is called behind) 
 launch_all.sh   : launch all testcases, but compilation should already be done
                   (run_test.sh is called behind) 

non_regression_test/testcases directory:

 tearing_limiter_199   : encompasses what is needed for tearing_limiter_199  test
 tearing_xpoint_302    : encompasses what is needed for tearing_xpoint_302   test
 tearing_xpoint_303    : encompasses what is needed for tearing_xpoint_3303   test
 get_testcase_data.sh  : get restart files from jorek.eu
 send_testcase_data.sh : upload restart files to jorek.eu

non_regression_test/job_scripts directory:

 helios          : directory with batch scripts for helios machine
 occigen	 : directory with batch scripts for occigen machine
 poincare        : directory with batch scripts for poincare machine



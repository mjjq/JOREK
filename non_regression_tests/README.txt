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

     git clone ssh://git.iter.org/stab/jorek.git jorek_git
     git checkout -b nrt feature/nrt

2) Download the reference test cases (restart files) from the
   jorek.eu web site 

     jorek_git/non_regression_tests/get_all_data.sh
  
3) Transfer the directory to Helios 
     
      rsync -av --progress jorek_git USER@helios.iferc-csc.org:
      ssh USER@helios.iferc-csc.org
      cd jorek_git

4) Load the needed module and put in your PATH the compiler you will use.

      source non_regression_tests/job_scripts/helios/env.sh

   Copy/paste a Makefile.inc that works for the machine you are using in the
   jorek directory.
 
      cp Make.inc/Makefile.helios_h5_linux_intel.inc Makefile.inc


5) Compile the jorek executables for the test "tearing_xpoint_303" 
   using the command:
      non_regression_tests/run_test.sh -p tearing_xpoint_303

4) Run the simulation using a batch script that calls:
$ non_regression_tests/run_test.sh -k -n tearing_limiter_199

   As an example, you can look at job_scripts/poincare/tearing199.job
   to design a batch script that does it on a parallel machine.

Optional) To send a new set of restart files to the web site that encompass
   the database :
$  cd non-regression-tests/testcase; ./send-testcase-data.sh tearing_limiter_199

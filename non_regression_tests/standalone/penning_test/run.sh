#!/bin/bash

# Exit on error
set -e

# Make sure it is compiled first (in subshell so we stay in this dir)
JOREK_DIR=../../..
(cd $JOREK_DIR && make penning_test)

# Run all tests in the subdir cases
for testcase in cases/*; do
   echo "-------------------------------------------------------------------------------"
   echo "- " $testcase
   echo "-------------------------------------------------------------------------------"
   # Run the test, and quit on the error code (due to set -e above)
   $JOREK_DIR/penning_test < $testcase
done

#!/bin/bash

# Exit on error
set -e

# Make sure it is compiled first (in subshell so we stay in this dir)
JOREK_DIR=../../..
(cd $JOREK_DIR && make penning_test)

# Run the test, and reuse the exit code
$JOREK_DIR/penning_test < penning_in

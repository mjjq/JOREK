#!/bin/bash

# This script download a jorek test case.
TESTNAME="$1"
if [ ! -f ${TESTNAME}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi

exit 0

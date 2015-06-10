#!/bin/bash

# This script upload a jorek test case.
TESTNAME="$1"
URL="http://jorek.eu/dav_nrt/"
#URL="http://localhost:8034/dav_nrt/"
if [ ! -f ${TESTNAME}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi
cd  ${TESTNAME}
curl -u nrt:nrt_21745XtL -T ${TESTNAME}.tgz  ${URL} 

exit 0

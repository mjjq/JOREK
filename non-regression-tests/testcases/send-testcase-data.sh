#!/bin/bash

# This script upload a jorek test case.
TESTDIR="$1"
TESTNAME=$(basename $TESTDIR)
URL="http://jorek.eu/dav_nrt/"
#URL="http://localhost:8034/dav_nrt/"
if [ ! -f ${TESTDIR}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi

echo "Creating tarball ${TESTDIR}.tgz"
cd  ${TESTDIR}
testcasedir=`readlink -f ${PWD}`
source ./settings.sh
pack_restart_files || exit 1

echo "Uploading ${TESTDIR}.tgz"
curl -s -u nrt:nrt_21745XtL -T ${TESTNAME}.tgz  ${URL} 
if [ $? -eq 0 ]; then
  printf "Success\n"
  exit 0
else
  printf "Failed\n"
  exit 1
fi


#!/bin/bash

# This script upload a jorek test case.
TESTDIR="$1"
URL="http://jorek.eu/dav_nrt/"
#URL="http://localhost:8054/dav_nrt/"
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
TESTNAME=$(basename $TESTDIR)
http_proxy=http://proxy:3128 curl -s -u nrt:nrt_21745XtL -T ${TESTNAME}.tgz  ${URL} 
if [ $? -eq 0 ]; then
  printf "Success\n"
  exit 0
else
  printf "Failed\n"
  exit 1
fi


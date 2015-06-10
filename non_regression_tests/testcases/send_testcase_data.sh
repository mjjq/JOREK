#!/bin/bash

# This script upload a jorek test case.
TESTDIR="$1"
URL="http://jorek.eu/dav_nrt/"
if [ ! -f ${TESTDIR}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi

echo "Creating tarball ${TESTDIR}.tgz"
cd  ${TESTDIR}
testcasedir=`readlink -f ${PWD}`
source ./settings.sh

testname=$(basename $testcasedir)
tar cvzf ${testname}.tgz begin.h5 end.h5                                           || exit 1

echo "Uploading ${TESTDIR}.tgz"
TESTNAME=$(basename $TESTDIR)
curl -s -u nrt:nrt_21745XtL -T ${TESTNAME}.tgz  ${URL}

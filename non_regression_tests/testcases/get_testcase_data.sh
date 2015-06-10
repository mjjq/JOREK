#!/bin/bash

# This script download a jorek test case.
TESTNAME="$1"
URL="http://jorek.eu/dav_nrt"
if [ ! -f ${TESTNAME}/settings.sh ]; then
  printf "This test name does not exist\n"
  exit 1
fi

cd  ${TESTNAME}
if [ ! -f ${TESTNAME}.tgz ]; then 
  echo "Downloading ${TESTNAME}.tgz"
  wget -q --user=nrt --password=nrt_21745XtL ${URL}/${TESTNAME}.tgz 
  returncode=$?
  if [ $returncode -ne 0 ]; then
    cat <<EOF
###################################################################################
  Failed to automatically download from web site.
  Please download reference data ${TESTNAME}.tgz yourself 
  from http:/jorek.eu/dav_nrt and
  copy it into testcases/${TESTNAME} directory.
  Launch this script again to decompress the archive.
###################################################################################
EOF
    exit 1
  fi
else 
    cat <<EOF
###################################################################################
 No downloading performed because ${TESTNAME}.tgz already
 exists, suppress it if you want to download:
   rm testcases/${TESTNAME}/${TESTNAME}.tgz
###################################################################################                             
EOF
    exit 1
fi

echo "Uncompress ${TESTNAME}.tgz"
tar xvzf ${TESTNAME}.tgz
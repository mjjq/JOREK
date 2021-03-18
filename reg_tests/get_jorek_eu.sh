#!/bin/bash

# This script download a jorek test environment and Makefile.inc
TESTNAME=`echo "$1" | sed -e 's|/$||'`
if [ -z "${DAV_URL}" ]; then
    DAV_URL="http://jorek.eu/dav_nrt"
fi

while [ $# -gt 0 ]; do 
wget -q --user=nrt --password=nrt_21745XtL ${DAV_URL}/$1 || exit 1
shift
done


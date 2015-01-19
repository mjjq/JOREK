#!/usr/bin/env bash

# -----------------------------------------------------------------
# script used to read BINARY restart file
# and write values in HDF5 restart file
# -----------------------------------------------------------------
if [ $# -lt 2 ]; then
   echo "==========================================================="
   echo " Syntaxe $0 RES_DIRECTORY iter_number "
   echo "    Example:  ./$0 DC_TEST_HDF5 89    "
   echo " Remark: rst_bin2hdf5 executable should be in RES_DIRECTORY"
   echo "==========================================================="
   echo " "
   exit
fi

# Input binary file
# -----------------
# RST_DIR="DC_TEST_HDF5_BIS_JENKINS199"
RST_DIR=$1

ITERNUMBER=$2
#                jorekxxx43.rst
RST_BIN=`printf "jorek%5.5d.rst" $d $ITERNUMBER`
printf "%5.5d \n" $d $ITERNUMBER > ./file.out

FILE=${RST_DIR}'/'${RST_BIN}

echo "==========================================="
echo " Reading BINARY restart file :  $FILE "
echo "==========================================="

if [ -f $FILE ]; then
    echo " FILE $FILE exist "
    cp $FILE $RST_BIN
else
    echo " File $FILE not exist "
    echo "  =>  so exit !!"
    exit 1
fi

# Environment
# -----------
echo " Running  " $RST_DIR/rst_bin2hdf5 $FILE
$RST_DIR/rst_bin2hdf5
rm ./file.out ./$RST_BIN





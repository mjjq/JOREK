#!/usr/bin/env bash

# -----------------------------------------------------------------
# script used to read HDF5 restart file
# and write values in BINARY restart file
# -----------------------------------------------------------------
if [ $# -lt 1 ]; then
   echo "==========================================================="
   echo " Syntaxe $0 RES_DIRECTORY iter_number "
   echo "    Example:  ./$0 DC_TEST_HDF5 89    "
   echo " Remark: rst_hdf52bin executable should be in RES_DIRECTORY"
   echo "==========================================================="
   echo " "
   exit
fi

# Input binary file
# -----------------
# RST_DIR="../DC_TEST_HDF5_BIS_JENKINS199/"
RST_DIR=$1

ITERNUMBER=$2
#                jorekxxxxx.h5
RST_h5=`printf "jorek%5.5d.h5" $d $ITERNUMBER`
printf "%5.5d \n" $d $ITERNUMBER > ./file.out

FILE=${RST_DIR}'/'${RST_h5}

echo "==========================================="
echo " Reading HDF5 restart file :  $FILE "
echo "==========================================="

if [ -f $FILE ]; then
    echo " FILE $FILE exist "
    cp $FILE $RST_h5
else
    echo " File $FILE not exist "
    echo "  =>  so exit !!"
    exit 1
fi

# Environment
# -----------
echo " Running  " $RST_DIR/rst_hdf52bin $FILE
$RST_DIR/rst_hdf52bin
rm -rf ./file.out ./$RST_h5




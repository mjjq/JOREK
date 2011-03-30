#!/bin/bash

#
# Purpose: Extract data from 'macroscopic_vars.dat' which is written by JOREK during the code run.
#
# Data: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

function usage() {
  echo ""
  echo "Usage: $0 <PREFIX> [<target_file>]"
  echo ""
  echo "  if <target_file> is omitted or '-', the output goes to STDOUT"
  echo ""
  echo "  e.g. extract_live_data ENERGIES energies.dat"
  echo "  e.g. extract_live_data N_TOR -"
  echo ""
}

# --- Evaluate command line parameters
if [ $# -lt 1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage
  exit
fi

PREFIX=$1
TARGET=$2

function extract() { grep "@${PREFIX}:" macroscopic_vars.dat | sed -e "s/^@${PREFIX}://"; }

if [ "$TARGET" == "-" ] || [ "$TARGET" == "" ]; then
  extract
else
  extract > ${TARGET}
fi

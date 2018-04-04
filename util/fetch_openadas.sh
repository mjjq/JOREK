#!/usr/bin/env bash
# Retrieve the OPEN-ADAS ADF11 datafiles for a specific symbol
#
# Author: Daan van Vugt <daanvanvugt@gmail.com>
# Date: 2018-03-28

set -o pipefail
set -o errtrace
set -o nounset
set -o errexit

declare -a ADF11_sets=(acd scd qcd xcd ccd plt prb pls prc)

function usage() {
  echo "Retrieve OPEN-ADAS ADF11 datafiles from the OPEN-ADAS website"
  echo "and store them in files (type_\$symbol.dat) in the current directory."
  echo "by default we try to download all possible sets in ADF11 with this suffix,"
  echo "which are $(tr '[:lower:]' '[:upper:]' <<< ${ADF11_sets[@]})."
  echo "See http://open.adas.ac.uk/adf11 for more information."
  echo ""
  echo "Usage: $0 suffix [suffix2] [...]"
  echo "  where suffix is an identifier such as 50_w. See the OPEN-ADAS website"
  echo "  for more information on the different datasets."
  echo ""
  echo "Examples:"
  echo "  $0 50_w"
  echo ""
  echo "Options:"
  echo "  -h, --help: show this message"
}

if [ "$#" -le 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  usage && exit
fi

for element in "$@"; do
  type=(${element//_/ }) # an array of 50 and w
  for set in "${ADF11_sets[@]}"; do
    if wget -q http://open.adas.ac.uk/download/adf11/$set$type/$set$element.dat -O $set$element.dat; then
      # check if it is a correct file
      if grep -q "You have an error" $set$element.dat; then
        rm $set$element.dat
      else
        echo "Downloaded $set$element.dat"
      fi
    fi
  done
done

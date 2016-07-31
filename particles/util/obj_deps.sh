#!/usr/bin/env bash
set -euf -o pipefail

# Check all of the .dep/*.d files for dependencies of $1
MAX_ITERATIONS=20 # mostly here to catch loops

# to_check contains names of .d files that still need to be checked
to_check=`echo "$1" | sed -e 's/^[.]obj/.dep/' -e 's/[.]o$/.d/'`
# deps contains a list of dependencies
deps=""

while [ -n "$to_check" ]; do
  to_check_new=""
  for depfile in $to_check; do
    stems=`grep -o ' .mod/[^.][^/ ]*[.]mod' "$depfile" |\
      sed -e 's/[.]mod//g' -e 's/[/ ]//g'`
    for stem in $stems; do
      if ! $(echo "$deps" | grep -qF ".obj/$stem.o"); then
        # if not already in dependency list
        # Add to dependency and check lists
        deps="$deps .obj/$stem.o"
        echo ".obj/$stem.o" # echo here because $deps is unset after the loop
        to_check_new="$to_check_new .dep/$stem.d"
      fi
    done
  done
  to_check="$to_check_new"
done

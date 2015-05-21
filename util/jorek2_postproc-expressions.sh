#!/bin/bash

# This is an auxilliary script which outputs some help information of jorek2_postproc
# such that it can be copy-pasted into the jorek.eu/wiki


out="./for-wiki.txt"

echo "^ Expression ^ Description ^" > $out

echo "expressions" | jorek2_postproc | grep -A 999 000001 | grep -B 999 "\-------" | head -n -1 | sed -e 's/^ *[0-9]* *//' | sed -e 's/$/|/' >> $out

echo "Output written to $out."

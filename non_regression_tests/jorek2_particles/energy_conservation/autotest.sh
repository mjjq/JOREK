#!/bin/bash
set -e

# How much we can be above the reference result
SAFETY_FACTOR=1.5 # this is large because of randomness in positions

# Run testcase
./run.sh -v jet_eq_autotest

# Look for output of autotest case
# Compare momentum
echo "Comparing momentum"
tail -n 1 jet_eq_autotest_out/momentum_results.txt | cut -d' ' -f1 | awk "{if(\$1 > $SAFETY_FACTOR*1.223953649275522189e-07) exit 1}" 
echo "Succes!"

# Compare energy
echo "Comparing energy"
tail -n 1 jet_eq_autotest_out/energy_results.txt | cut -d' ' -f1 | awk "{if(\$1 > $SAFETY_FACTOR*4.059586000693115845e-14) exit 1}" 

echo "Succes!"

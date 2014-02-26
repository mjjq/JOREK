#!/bin/bash

out="./for-wiki.txt"

echo "^ Expression ^ Description ^" > $out

echo "expressions" | jorek2_postproc | grep -A 999 000001 | grep -B 999 "\-------" | head -n -1 | sed -e 's/^ *[0-9]* *//' | sed -e 's/$/|/' >> $out

echo "Output written to $out."

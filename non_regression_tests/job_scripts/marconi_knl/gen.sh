#!/bin/bash
LIST="kink_circ_710.job \
tear_circ_199.job	\
tear_circ_303_fft.job	\
tear_circ_303_gears.job	\
tear_circ_303.job	\
tear_circ_303_neo.job	\
tear_circ_303_tauic.job	\
tear_circ_303_tgnum.job \
tear_circ_333.job"

for i in $LIST; do NAME=${i/.job/}; sed  "s/_CASE_/$NAME/g" case.skel > $i; done
chmod +x $LIST

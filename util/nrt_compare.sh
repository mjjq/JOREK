# Location of 'util' directory 
UTILDIR=$(readlink -f `dirname $0`)
EXTSCRIPT=$UTILDIR/extract_live_data.sh
SED=sed
RPTFILE=report.txt

if [ $# -lt 3 ] || [ ! -d "$2" ] || [ ! -d "$3" ]; then
  printf "Usage: $0 <model_no> <Directory 1> <Directory 2>\n"
  exit 0
fi

if [ "$1" = "199" ]; then 
  THR=1e-2
  BEGL=23
  ENDL=75
  COLS="1,3,5"
fi
if [ "$1" = "302" ]; then 
  THR=1e-2
  BEGL=150
  ENDL=290
  COLS="1,3,5"
fi
if [ "$1" = "303" ]; then 
  THR=1e-2
  BEGL=350
  ENDL=900
  COLS="1,3,5"
fi
(cd $2; $EXTSCRIPT growth_rates) > o1
(cd $3; $EXTSCRIPT growth_rates) > o2
NBL1=`wc -l < o1`
NBL2=`wc -l < o2`
if [ $NBL1 -gt $NBL2 ]; then ((NBL=NBL2)); else ((NBL=NBL1)); fi
if [ $NBL  -lt $ENDL ]; then ((ENDL=NBL)); fi
for i in $(seq 1 2); do
  ${SED} -n "1p;${BEGL},${ENDL}p" < o$i | sed "s/  */ /g" | cut -d" " -f"${COLS}" > f$i
done
numdiff f1 f2 -r $THR 1> $RPTFILE
RET=$?
TOTLCOMP=`wc -l < f1`
if [ $TOTLCOMP -lt 3 ]; then
    printf "FAILED (Comparison can not be done, simulation not advanced enough)\n"
   exit 1
fi
if [ $RET -eq 0 ]; then
  printf "OK  (nb lines compared: $((TOTLCOMP-1))/$((ENDL-BEGL+1)))\n"
else
  printf "FAILED (lines compared : $TOTLCOMP)\n"
  cat $RPTFILE
fi
printf "\n# gnuplot commands to look at growth rates that have been compared :\n"
printf "    set key autotitle columnhead;\n"
printf "    set auto; plot 'f1' u 1:2 ls 1, 'f2' u 1:2 ls 4; pause -1\n"
printf "    set auto; plot 'f1' u 1:3 ls 3, 'f2' u 1:3 ls 6\n"
#rm -f o1 o2 
exit $RET
# Gnuplot plotting script used by plot_live_data script

#
# Purpose: Plot energies and growth rates during or after a JOREK code run
#
# Data: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

ps=<ps>
qtty="<qtty>"
ncols=<ncols>
logy=<logy>

if ( ps==1 ) set term postscript enhanced color
if ( ps==1 ) set output qtty.'.ps'
if ( logy==1 ) set log y
set key outside
set title qtty
plot for [i=2:ncols+1] qtty.'.dat' u 1:i w lp t columnhead(i)
pause mouse key

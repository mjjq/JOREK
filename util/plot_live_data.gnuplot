# Gnuplot plotting script used by plot_live_data script

#
# Purpose: Plot energies and growth rates during or after a JOREK code run
#
# Date: 2011-03-30
# Author: Matthias Hoelzl, IPP Garching
#

ps=<ps>
qtty="<qtty>"
ncols=<ncols>
logy=<logy>
xlabel="<xlabel>"
ylabel="<ylabel>"

if ( ps==1 ) set term postscript enhanced color
if ( ps==1 ) set output qtty.'.ps'
if ( logy==1 ) set log y
set key outside
set title qtty
set xlabel xlabel
set ylabel ylabel
plot for [i=2:ncols+1] qtty.'.dat' u 1:i w lp t columnhead(i)

if ( ps==0 ) print ''
if ( ps==0 ) print 'Hints:'
if ( ps==0 ) print '* Use right mouse button to zoom'
if ( ps==0 ) print '* Press "U" to unzoom again'
if ( ps==0 ) print '* Press "E" to update the data'
if ( ps==0 ) print '* Click with left mouse button to exit'
if ( ps==0 ) print ''
if ( ps==0 ) pause mouse button1

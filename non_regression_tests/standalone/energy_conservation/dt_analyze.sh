#!/bin/bash

# Gather results from different folders
rm -f energy_dt momentum_dt
for dir in *_t*[0-9]_out; do
   dt=`echo $dir | sed 's/[^0-9.]*//g'`
   echo -n "$dt " >> energy_dt
   echo -n "$dt " >> momentum_dt
   tail -n 1 $dir/energy_results.txt >> energy_dt
   tail -n 1 $dir/momentum_results.txt >> momentum_dt
done

# TODO fit automatically
# Plot results
cat <<EOF | gnuplot
   set logscale xy
   set xlabel "dt"
   set terminal png
   set style data points

   set ylabel "Growth rate abs((end-begin)/begin-1)"
   set out 'energy_growth.png'
   plot 'energy_dt' u 1:2 t "average", 'energy_dt' u 1:3 t "stddev", 'energy_dt' u 1:4 t "minimum", 'energy_dt' u 1:5 t "maximum", 10**(-13.0)/sqrt(x) t 'dt^(-0.5)'
   set out 'momentum_growth.png'
   set key bottom right
   plot 'momentum_dt' u 1:2 t "average", 'momentum_dt' u 1:3 t "stddev", 'momentum_dt' u 1:4 t "minimum", 'momentum_dt' u 1:5 t "maximum", 0.01*x*x t 'dt^2'
   set key top right

   set ylabel "Relative energy change"
   set out 'energy_variance.png'
   plot 'energy_dt' u 1:6 t "normalized standard deviation", 10**(-15.5)/sqrt(x) t 'dt^(-0.5)'
   set out 'momentum_variance.png'
   plot 'momentum_dt' u 1:6 t "normalized standard deviation", 10**(-1.5)*x*x t 'dt^2'

   set ylabel "1 - Min/Max ratio"
   set out 'energy_minmax.png'
   plot 'energy_dt' u 1:7 t "highest"
   set out 'momentum_minmax.png'
   set key bottom right
   plot 'momentum_dt' u 1:7 t "highest", 10**(-0.5)*x*x t 'dt^2'
   set key top right
EOF

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
   set terminal cairolatex size 10cm,7.5cm
   set style data points

   set ylabel "Growth rate abs((end-begin)/begin-1)"
   set out 'energy_growth.png'
   plot 'energy_dt' u 1:2 t "average", 'energy_dt' u 1:3 t "stddev", 'energy_dt' u 1:4 t "minimum", 'energy_dt' u 1:5 t "maximum", 10**(-13.0)/sqrt(x) t 'dt^(-0.5)'
   set out 'momentum_growth.png'
   set key bottom right
   plot 'momentum_dt' u 1:2 t "average", 'momentum_dt' u 1:3 t "stddev", 'momentum_dt' u 1:4 t "minimum", 'momentum_dt' u 1:5 t "maximum", 0.01*x*x t 'dt^2'
   set key top right

   set format y "\$10^{%L}\$"
   set ylabel "Energy stddev"
   set out 'energy_variance.tex'
   fit a/sqrt(x) 'energy_dt' u 1:6 via a
   plot 'energy_dt' u 1:6 t "normalized standard deviation", a/sqrt(x) t 'fit \$\\Delta t^{-0.5}\$'
   set ylabel "Momentum stddev"
   set out 'momentum_variance.tex'
   fit a*x**2 'momentum_dt' u 1:6 via a
   set key top left
   plot 'momentum_dt' u 1:6 t "normalized standard deviation", a*x*x t 'fit \$\\Delta t^2\$'

   set ylabel "1 - Min/Max ratio"
   set out 'energy_minmax.png'
   plot 'energy_dt' u 1:7 t "highest"
   set out 'momentum_minmax.png'
   set key bottom right
   plot 'momentum_dt' u 1:7 t "highest", 10**(-0.5)*x*x t 'dt^2'
   set key top right
EOF

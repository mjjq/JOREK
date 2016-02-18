set logscale xy
set xlabel '$\Delta t$ [JOREK time]'
set ylabel '$|x-x_{analytical}|$ [m]'
set term cairolatex size 10cm,7.5cm
name='dt_results/11_multi_element_polar'
set out name.'.tex'
fit a*x**2 name.'.dat' u 5:7 via a
set key bottom right
p name.'.dat' u 5:7 w p t 'position error', a*x**2 w l t sprintf('%g $\Delta t^2$', a)

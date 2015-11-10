name='n_radial_results/n_radial'
set logscale xy
set xlabel 'N_{radial} and N_{poloidal}'
set ylabel '|x-x_{analytical}| [m]'
set term pngcairo
set out name.'.png'
fit a*x**-4 name.'.dat' u 3:7 via a
set key bottom right
p name.'.dat' u 3:7 w p t 'position error', a*x**-4 w l t sprintf('%g N^{-4}', a)

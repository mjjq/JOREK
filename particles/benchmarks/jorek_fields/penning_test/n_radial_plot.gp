set logscale xy
set xlabel '$N_{radial}$ and $N_{poloidal}$'
set ylabel '$|x-x_{analytical}|$ [m]'
set term cairolatex size 10cm,7.5cm
set out 'n_radial_results/n_radial.tex'
fit a*x**-4 'n_radial_results/n_radial.dat' u 3:7 via a
set key top right
set xr [10:200]
p 'n_radial_results/n_radial.dat' u 3:7 w p t 'position error', a*x**-4 w l t 'fit' #sprintf('%g $N^{-4}$', a)

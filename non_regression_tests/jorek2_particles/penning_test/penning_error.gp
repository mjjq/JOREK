set logscale xy
set xlabel '{/Symbol D}t [JOREK time]'
set ylabel '|x-x_{analytical}| [m]'
#set term pngcairo
#set out 'penning_error.png'
fit a*x**2 'penning_error.out' u 1:2 via a
set key bottom right
p 'penning_error.out' u 1:2 w p t 'position error', a*x**2 w l t sprintf('%g {/Symbol D}t^2', a)

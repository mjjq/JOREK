set logscale xy
set xlabel '{/Symbol D}t [JOREK time]'
set ylabel '|x-x_{analytical}| [m]'
set term pngcairo
set out name.'.png'
fit a*x**2 name.'.dat' u 5:7 via a
set key bottom right
p name.'.dat' u 5:7 w p t 'position error', a*x**2 w l t sprintf('%g {/Symbol D}t^2', a), name.'.dat' u 5:9 w l t 'reference'

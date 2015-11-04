set xlabel "R"
set ylabel "Z"
p 'positions_RZ.dat' every ::1 u 2:3 t 'Particle #1' w l, \
  'gc_RZ.dat' every ::1 u 2:3 t 'Particle #1 GC' w l
  #'positions_RZ.dat' every ::1 u 4:5 t 'Particle #2' w p, \
  #'positions_RZ.dat' every ::1 u 6:7 t 'Particle #3' w p, \

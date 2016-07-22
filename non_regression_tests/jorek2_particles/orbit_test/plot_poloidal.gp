set xlabel "R"
set ylabel "Z"
p 'gc_RZ.dat' every ::1 u 2:3 t 'Particle #1 GC' w l, \
  'gc_RZ.dat' every ::1 u 5:6 t 'Particle #2 GC' w l, \
  'gc_RZ.dat' every ::1 u 8:9 t 'Particle #3 GC' w l

#p 'positions_RZ.dat' every ::1 u 2:3 t 'Particle #1' w l, \
  #'positions_RZ.dat' every ::1 u 4:5 t 'Particle #2' w l, \
  #'positions_RZ.dat' every ::1 u 6:7 t 'Particle #3' w l, \

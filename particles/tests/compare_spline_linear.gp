set multiplot layout 2,1
set xr [0:74]
plot 'mc_time.txt' using 2:1:3 with image t ''
plot 'coronal_time.txt' using 2:1:3 with image t ''

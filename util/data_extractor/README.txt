This code can be used to extract data on a rectangular grid from a JOREK .h5 file.
The code does not need to be compiled with the rest of JOREK, it is independent.
The RZ-rectangular grid is defined up to the maxima of the JOREK domain.
If points are outside the JOREK domain, all variables are set to zero, expect psi_norm, which is set to 100.0

This code should work with any recent hdf5 JOREK file, but will not work with old hdf5 formats that are missing some equilibrium information like psi_axis and psi_bnd and central_density.
It needs to be compiled with an mpi compiler. The Makefile as it is here should work on Marconi.

Then, run the code to extract the data.
Data can be extracted in .txt of .h5 format.
Many variables are available, including derived variables like the magnetic and the electric field.
An option is also available to create an image of the data.

The code is parallelised, and the best way to use it, is first to create the "pixel" file.
This saves the grid coords from JOREK that you want. After this is done, you can run for any JOREK file (from the same simulation run), for any phi angle, using this "pixel file", which makes it much much faster.
And of course you can run in parallel, so "mpirun -np 16 ./etc"

Once compiled, here is an example using 16 mpi processes, and a grid resolution of 0.005m, with a JOREK restart file called jorek01234.h5:

# This will give you a list of the options
./fortran_process -h

# This saves the grid data
mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -save_pixels -resolution 0.005

# This saves the data, note that you can include many variables in a single command (option -h lists all the available variables)
# First we do this on 4 phi angles for file jorek01234.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.00000 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.31415 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.62830 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.94245 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5

# And you can repeat this for any another file, say jorek01235.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01235.h5 -use_pixel_file -phi 0.00000 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01235.h5 -use_pixel_file -phi 0.31415 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01235.h5 -use_pixel_file -phi 0.62830 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5
mpirun -np 16 ./fortran_process -jorek_file jorek01235.h5 -use_pixel_file -phi 0.94245 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data.h5




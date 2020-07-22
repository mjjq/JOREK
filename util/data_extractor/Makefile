# --- Makefile for jorek processing with hdf5

# --- Compiler
#FCOMP = pgf90
#FCOMP = f95
FCOMP = mpiifort
#FCOMP = mpif90

# --- HDF5 libs
#HDF5_HOME = /usr/local/depot/hdf5-1.8.5
#HDF5_HOME = /usr/local/depot/hdf5-1.8.19-ifort12
#HDF5_INC  = $(HDF5_HOME)/include
#HDF5_LIB  = $(HDF5_HOME)/lib/libhdf5_fortran.a -L$(HDF5_HOME)/lib -lhdf5
#HDF5_HOME = $(HDF5_HOME)
HDF5_INC  = $(HDF5_HOME)/include
HDF5_LIB  = $(HDF5_HOME)/lib/libhdf5_fortran.a -L$(HDF5_HOME)/lib -lhdf5

all:
	$(FCOMP) -I$(HDF5_INC) -c hdf5_io_module.f90
	$(FCOMP) -I$(HDF5_INC) -c fortran_process_jorek.f90
	$(FCOMP) -I$(HDF5_INC) fortran_process_jorek.o hdf5_io_module.o -o fortran_process $(HDF5_LIB)
	rm -f *.o *.mod


image:
	$(FCOMP) create_image.f90 -o create_image
	rm -f *.o 
	

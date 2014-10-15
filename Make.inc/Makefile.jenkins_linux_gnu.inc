# Confiduration file for jorek
# model directory
MODEL = model199

# Fortran compiler 
FC     = mpif90 -f90=gfortran
CC     = mpicc

FFLAGS_OMP       = -fopenmp 
FFLAGS          := $(FFLAGS) -malign-double -g  
FFLAGS          := $(FFLAGS) -fdefault-real-8 -fdefault-double-8 
FFLAGS          := $(FFLAGS) -ffree-line-length-none
FFLAGS          := $(FFLAGS) -x f95-cpp-input
## TODO find equivalents of this ifort options
#FFLAGS          := $(FFLAGS)  -xT -align 
#FFLAGS          := $(FFLAGS) -ffpe-trap=invalid,zero,overflow,underflow -fbounds-check 

FFLAGS_FIXEDFORM = $(FFLAGS) -ffixed-form  -fbounds-check $(FFLAGS_OMP) 
FFLAGS_NOBOUNDS  = $(FFLAGS) -ffixed-form                 $(FFLAGS_OMP) 
FFLAGS          := $(FFLAGS)               $(FFLAGS_OMP) 

#FFLAGS_NO_OMP   = -fpp -r8 -i4 -O2 -xT -align
FFLAGS          := $(FFLAGS) #-DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0
CFLAGS          := -DSTAN_FLAG 

#NOHDF5FLAG  = -DNOHDF5

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 0
USE_MUMPS = 0
USE_WSMP   = 0
USE_HDF5  = 1

#Scotch library
SCOTCH_HOME  = /builds/jorek_tools/scotch
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# PASTIX
PASTIX_HOME      = /builds/jorek_tools/pastix/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`

#HDF5
HDF5_HOME   = /usr
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a \
	      $(HDF5_HOME)/lib/libhdf5_hl.a \
	      $(HDF5_HOME)/lib/libhdf5_fortran.a \
	      $(HDF5_HOME)/lib/libhdf5.a -lz \

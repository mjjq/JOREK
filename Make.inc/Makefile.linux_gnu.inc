# Configuration file for jorek
# model directory
MODEL = model303

# Fortran compiler 
FC     = mpifort
CC     = mpicc 
CXX    = mpicxx 

DEBUG = 1

FFLAGS          := $(FFLAGS) -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0 -DFUNNELED -DWORLDWAR2 -fdefault-real-8 -fdefault-double-8
FLAGS += -O3 -msse2 -march=native

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 1
USE_MUMPS = 1
USE_WSMP   = 0
USE_FFTW = 0
USE_HDF5 = 1

#Scotch library
SCOTCH_HOME  = $(HOME)/libs/scotch_5.1.12_esmumps
LIB_SCOTCH   = -L$(SCOTCH_HOME)/libscotchs -lscotch -lscotcherr $(HOME)/libs/scotch_5.1.12_esmumps/src/esmumps/libesmumps.a
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME = $(HOME)/libs/MUMPS
LIB_MUMPS  = $(MUMPS_HOME)/lib/libdmumps.a $(MUMPS_HOME)/lib/libmumps_common.a $(LIB_SCOTCH)
LIB_ZMUMPS = $(MUMPS_HOME)/lib/libzmumps.a 
INC_MUMPS  = $(MUMPS_HOME)/include/
ORDLIB     = $(MUMPS_HOME)/lib/libpord.a

# PASTIX
PASTIX_HOME      = $(HOME)/libs/pastix4492/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`

#HIPS
LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#WSMP
LIB_WSMP = -L$(WSMP_HOME) -lpwsmp64

# LIBFFTW
#LIBFFTW   = -L/usr/lib/x86_64-linux-gnu/ -lfftw3
#INCFFTW   = -L/usr/include/

# HDF5
# Install using apt install hdf5-tools libhdf5-openmpi-dev libhdf5-serial-dev
HDF5LIB     = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
HDF5INCLUDE = /usr/include/hdf5/openmpi

LIBLAPACK = /usr/lib/liblapack.a -lblas

SCALAP    = /usr/lib/libscalapack-openmpi.a #-lscalapack



# BLACS for ubuntu scalapack
BLACS_HOME = /usr/lib
BLACS       = $(BLACS_HOME)/libblacsCinit-openmpi.a     \
              $(BLACS_HOME)/libblacsF77init-openmpi.a   \
              $(BLACS_HOME)/libblacs-openmpi.a          \
              $(BLACS_HOME)/libblacsF77init-openmpi.a   \
              $(BLACS_HOME)/libblacsCinit-openmpi.a     \
              $(BLACS_HOME)/libblacs-openmpi.a

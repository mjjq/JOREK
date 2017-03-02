# Confiduration file for jorek
# model directory
MODEL = model199

# Fortran compiler 
FC     = mpif90
CC     = mpicc

FFLAGS_OMP       = -fopenmp 
FFLAGS           = -g -O0 
FFLAGS          := $(FFLAGS) -fdefault-real-8
FFLAGS          := $(FFLAGS) -ffree-line-length-none
FFLAGS          := $(FFLAGS) -x f95-cpp-input
## TODO find equivalents of this ifort options
#FFLAGS          := $(FFLAGS) -traceback -xT -align 

FFLAGS_FIXEDFORM = $(FFLAGS) -ffixed-form  -fbounds-check $(FFLAGS_OMP) 
FFLAGS_NOBOUNDS  = $(FFLAGS) -ffixed-form                 $(FFLAGS_OMP) 
FFLAGS          := $(FFLAGS)               -fbounds-check $(FFLAGS_OMP) 

#FFLAGS_NO_OMP   = -fpp -O2 -xT -align
FFLAGS          := $(FFLAGS) -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 0
USE_MUMPS  = 0
USE_WSMP   = 0
USE_FFTW   = 0
USE_HDF5   = 1

#Scotch library
SCOTCH_HOME  = /home/ITER/vanvugd/source/scotch_5.1.12
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME ?= 
LIB_MUMPS  = $(MUMPS_HOME)/lib/libdmumps.a
LIB_ZMUMPS = $(MUMPS_HOME)/lib/libzmumps.a
INC_MUMPS  = $(MUMPS_HOME)/include/
ORDLIB     = $(MUMPS_HOME)/lib/libpord.a

# PASTIX
PASTIX_HOME      = /home/ITER/vanvugd/source/pastix_5.2.2.22/install
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
LIBFFTW   = -L$(FFTWDIR)/lib -lfftw3
INC_FFTW  = -I$(FFTWDIR)/include

#LAPACK
#LIBLAPACK = -L$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#SCALAP    = -L$(MKL) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 

HDF5_HOME=$(H5DIR)
HDF5INCLUDE=$(HDF5_HOME)/include
HDF5LIB=$(HDF5_HOME)/lib/libhdf5hl_fortran.a \
        $(HDF5_HOME)/lib/libhdf5_hl.a        \
        $(HDF5_HOME)/lib/libhdf5_fortran.a   \
        $(HDF5_HOME)/lib/libhdf5.a           \
	$(EBROOTZLIB)/lib/libz.a


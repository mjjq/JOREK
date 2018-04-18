# Configuration file for jorek

# model directory
MODEL = model303

# Fortran compiler 
FC               = mpiifort -mt_mpi
CC               = mpiicc   -mt_mpi
CXX              = mpiicpc  -mt_mpi
FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4 -DUSE_BLOCK -DFUNNELED -DWORLDWAR2 #-DCOMPARE_ELEMENT_MATRIX
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -xHost -align -O2

# --- DEBUGGING
#   --- Debug symbols for debuggers
DEBUGFLAGS  = #-g -p
#DEBUGFLAGS += -warn all,nounused -check all,noarg_temp_created -debug all -debug-parameters -fstack-security-check -traceback

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

FFLAGS_NO_OMP    = $(FFLAGS) $(DEBUGFLAGS)

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 0
USE_MUMPS  = 0
USE_WSMP   = 0
USE_FFTW = 1
USE_HDF5 = 1

# PASTIX
LIB_PASTIX       = `$(PASTIX_HOME)/install/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/install/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/install/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/install/pastix-conf --incs`

# LIBFFTW
LIBFFTW   = $(FFTW_HOME)/lib/libfftw3.a
INC_FFTW = -I$(FFTW_HOME)/include

#HDF5
#HDF5_HOME   = ${HDF5}
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a \
              $(HDF5_HOME)/lib/libhdf5_hl.a \
              $(HDF5_HOME)/lib/libhdf5_fortran.a \
              $(HDF5_HOME)/lib/libhdf5.a -lz

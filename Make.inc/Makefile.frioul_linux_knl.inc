# Makefile.inc file for compiling JOREK on the Marconi-Fusion Supercomputer
#
# Make sure to set environment variables and load modules as described in 
#   http://jorek.eu/wiki/doku.php?id=marconi-fusion_compile_knl_intelmpi_i17&do=

# --- Select physics model
MODEL             = model303

# --- Compiler and options
FC                = mpif90 -mt_mpi
CC                = mpicc -mt_mpi
CXX               = mpicxx -mt_mpi
FFLAGS            = -i4 -DUSE_BLOCK #-DFUNNELED #RECENT_PASTIX#-DWORLDWAR2
FFLAGS           += -xmic-avx512 -O2 -fp-model fast=2

# --- Various switches
USE_PASTIX        = 1
USE_PASTIX_MURGE  = 0
USE_HDF5          = 1
USE_FFTW          = 1
USE_MUMPS         = 0

# --- PaStiX solver
PASTIX_LIB ?= "<<<<path-to-pastix-library-directory>>>>"
LIB_PASTIX        = `$(PASTIX_LIB)/pastix-conf --libs` -L${MKL_LIB}
LIB_PASTIX_MURGE  = `$(PASTIX_LIB)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS   = `$(PASTIX_LIB)/pastix-conf --blas`
INC_PASTIX        = `$(PASTIX_LIB)/pastix-conf --incs`

# --- FFTW Library
  #Through INTEL MKL 
LIBFFTW           = -L$(MKL_LIB_DIR) -lfftw3xf_intel
INC_FFTW          = -I$(MKL_INC_DIR)/fftw

# --- Usually not to be changed
FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NO_OMP     = $(FFLAGS) $(DEBUGFLAGS)

HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

# --- MUMPS solver
MUMPS_HOME ?= "<<<<path-to-your-Mumps-directory>>>>>"
MUMPS_LIB ?= $(MUMPS_HOME)/lib
MUMPS_INC ?= $(MUMPS_HOME)/include
LIB_MUMPS = -L$(MUMPS_LIB) -ldmumps -lmumps_common  -lpord #-lgfortran
INC_MUMPS = $(MUMPS_INC)                                                        
SCALAP    = -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

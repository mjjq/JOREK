# Makefile.inc file for compiling JOREK on the Marconi-Fusion Supercomputer
#
# Make sure to set environment variables and load modules as described in 
#   http://jorek.eu/wiki/doku.php?id=marconi-fusion_compile_knl_intelmpi_i17&do=

# --- Select physics model
MODEL             = model303

# --- Compiler and options
FC                = mpiifort -mt_mpi
CC                = mpiicc   -mt_mpi
CXX               = mpiicpc  -mt_mpi
FFLAGS            = -i4 -DUSE_BLOCK -DFUNNELED #RECENT_PASTIX#-DWORLDWAR2
FFLAGS           += -axmic-avx512 -qopt-prefetch=5 -O2

# --- Various switches
USE_PASTIX        = 1
USE_PASTIX_MURGE  = 0
USE_HDF5          = 1
USE_FFTW          = 1
USE_MUMPS         = 1

# --- PaStiX solver
PASTIX_LIB ?= "<<<<path-to-pastix-library-directory>>>>"
LIB_PASTIX        = `$(PASTIX_LIB)/pastix-conf --libs` 
LIB_PASTIX_MURGE  = `$(PASTIX_LIB)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS   = `$(PASTIX_LIB)/pastix-conf --blas`
INC_PASTIX        = `$(PASTIX_LIB)/pastix-conf --incs`

# --- FFTW Library
LIBFFTW           = $(FFTW_LIB)/libfftw3.a
INC_FFTW          = -I$(FFTW_INC)

# --- Usually not to be changed
FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NO_OMP     = $(FFLAGS) $(DEBUGFLAGS)

HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L$(SZIP_LIB) -lsz -lz

# --- MUMPS solver
MUMPS_HOME ?= "<<<<path-to-your-Mumps-directory>>>>>"
MUMPS_LIB ?= $(MUMPS_HOME)/lib
MUMPS_INC ?= $(MUMPS_HOME)/include
LIB_MUMPS = -L$(MUMPS_LIB) -ldmumps -lmumps_common  -lpord #-lgfortran
INC_MUMPS = $(MUMPS_INC)                                                        
SCALAP    = -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

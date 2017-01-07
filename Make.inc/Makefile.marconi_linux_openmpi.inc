# Makefile.inc file for compiling JOREK on the Marconi-Fusion Supercomputer
#
# * Make sure to set PASTIX_HOME below to the correct directory.
# * Make sure that you have the following modules loaded: intel intelmpi mkl fftw hdf5 szip

# --- Select physics model
MODEL             = model199

# --- Compiler and options
FC                = mpif90 #INTELMPI#-mt_mpi
CC                = mpicc  #INTELMPI#-mt_mpi
FFLAGS_OMP        = -qopenmp 
FFLAGS            = -r8 -i4 -DUSE_BLOCK -DFUNNELED #RECENT_PASTIX#-DWORLDWAR2
FFLAGS           += -fpp
FFLAGS           += -xHost -O2

# --- Debugging options
#DEBUGFLAGS        = -g
#DEBUGFLAGS      += -warn all,nounused
#DEBUGFLAGS      += -check all,noarg_temp_created
#DEBUGFLAGS      += -debug all
#DEBUGFLAGS      += -debug-parameters
#DEBUGFLAGS      += -fstack-security-check
#DEBUGFLAGS      += -ftrapuv

# --- Various switches
USE_PASTIX        = 1
USE_PASTIX_MURGE  = 0
USE_HDF5          = 1
USE_FFTW          = 1
USE_MUMPS         = 1

# --- PaStiX solver
PASTIX_LIBDIR ?= "<<<<path-to-pastix-library-directory>>>>"
LIB_PASTIX        = `$(PASTIX_LIBDIR)/pastix-conf --libs` -L$(MKL_LIB)
LIB_PASTIX_MURGE  = `$(PASTIX_LIBDIR)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS   = `$(PASTIX_LIBDIR)/pastix-conf --blas`
INC_PASTIX        = `$(PASTIX_LIBDIR)/pastix-conf --incs`

# --- FFTW Library
LIBFFTW           = /lib64/libfftw3.a
INC_FFTW          = -I/usr/include

# --- Usually not to be changed
FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NO_OMP     = $(FFLAGS) $(DEBUGFLAGS)

HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a \
              $(HDF5_HOME)/lib/libhdf5_hl.a \
              $(HDF5_HOME)/lib/libhdf5_fortran.a \
              $(HDF5_HOME)/lib/libhdf5.a \
              $(SZIP_LIB)/libsz.a -lz

# --- MUMPS solver
MUMPS_HOME ?= "<<<<path-to-your-Mumps-directory>>>>>"
MUMPS_LIB ?= $(MUMPS_HOME)/lib
MUMPS_INC ?= $(MUMPS_HOME)/include
LIB_MUMPS = -L$(MUMPS_LIB) -ldmumps -lmumps_common  -lpord #-lgfortran
INC_MUMPS = $(MUMPS_INC)                                                        
SCALAP    = -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

# Configuration file for jorek

# model directory
MODEL = model302

# Fortran compiler 
FC     = mpif90 -f90=ifort
CC     = mpicc
CXX    = mpicxx -cxx=icpc

FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -traceback -xSSSE3 -align
FFLAGS          := $(FFLAGS) -O2

# --- DEBUGGING
#   --- Debug symbols for debuggers
DEBUGFLAGS  = -g
#   --- Check array boundaries
DEBUGFLAGS += -check bounds
#   --- Generate and check interfaces (__genmod files)
#DEBUGFLAGS += -gen-interfaces -warn interfaces
#   --- Check many other things
#DEBUGFLAGS += -warn all,nounused -check all,noarg_temp_created -debug all -debug-parameters -fstack-security-check -ftrapuv
#   --- Even more checks and warnings
#DEBUGFLAGS += -warn all -check all

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

#FFLAGS_NO_OMP    = -fpp -r8 -i4 -O2 -xSSSE3 -align
FFLAGS          := $(FFLAGS) -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 1
USE_MUMPS = 0
USE_WSMP   = 0

# Directory containing the user-compiled libraries (e.g., MUMPS, PASTIX, ...)
LIBDIR = $(HOME)/mylib/

#Scotch library
#SCOTCH_HOME ?= $(HOME)/scotch_5.1
#LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
#INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME = $(LIBDIR)/MUMPS_4.9.2
LIB_MUMPS  = -L$(MUMPS_HOME)/lib/ -ldmumps -lmumps_common
LIB_ZMUMPS = $(MUMPS_HOME)/lib/libzmumps.a
INC_MUMPS  = $(MUMPS_HOME)/include/
ORDLIB    = $(MUMPS_HOME)/lib/libpord.a #$(LIBDIR)/metis-4.0/libmetis.a $(LSCOTCH)


# PASTIX
PASTIX_HOME     ?= $(LIBDIR)/pastix_release_3131/install/
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`


#HIPS
LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#WSMP
LIB_WSMP = -L$(WSMP_HOME) -lpwsmp64

#FFTW
LIBFFTW = $(FFTW_ROOT)/lib/libfftw3.a

#LAPACK
LIBLAPACK = -L$(MKLPATH) -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5

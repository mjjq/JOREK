# Configuration file for jorek

# model directory
MODEL = model199

# Fortran compiler 
FC                = mpiifort -mt_mpi
CC                = mpiicc   -mt_mpi
CXX               = mpiicpc  -mt_mpi

FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -traceback -xSSSE3 -align
FFLAGS          := $(FFLAGS) -O0     # Optimization level

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

#Scotch library
SCOTCH_HOME ?= $(HOME)/scotch_5.1
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME = $(HOME)/intel/mumps/MUMPS_4.7.3
LIB_MUMPS  = $(MUMPS_HOME)/lib/libdmumps.a
LIB_ZMUMPS = $(MUMPS_HOME)/lib/libzmumps.a
INC_MUMPS  = $(MUMPS_HOME)/include/
ORDLIB     = $(MUMPS_HOME)/lib/libpord.a

# PASTIX
PASTIX_HOME     ?= $(HOME)/pastix
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
#LIBFFTW   = $(HOME)/intel/fftw/install/lib/libfftw3.a
LIBFFTW   = -L/opt/cluster/fftw/3.2.2/intel-11.1-069/intel-4.0.0.025 -lfftw3

#LIBLAPACK = $(HOME)/intel/lapack/LAPACK/lapack_LINUX.a -lblas
LIBLAPACK = -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core 
SCALAP    = -L$(HOME)/intel/scalapack/SCALAPACK -lscalapack
#SCALAP    = -lmkl_scalapack_lp64
BLACS_HOME = $(HOME)/intel/blacs/BLACS/LIB
BLACS      = $(BLACS_HOME)/blacs_MPI-LINUX-0.a		\
	     $(BLACS_HOME)/blacsF77init_MPI-LINUX-0.a 	\
	     $(BLACS_HOME)/blacsCinit_MPI-LINUX-0.a   	\
	     $(BLACS_HOME)/blacs_MPI-LINUX-0.a 
#BLACS     = -lmkl_blacs_intelmpi20_lp64

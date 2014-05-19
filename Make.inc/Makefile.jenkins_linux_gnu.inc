# Confiduration file for jorek
# model directory
MODEL = model199

# Fortran compiler 
FC     = mpif90 -f90=gfortran
CC     = mpicc

FFLAGS_OMP       = -fopenmp 
FFLAGS           = -O0 # -g  -ffpe-trap=invalid,zero,overflow -traceback -fbounds-check 
FFLAGS          := $(FFLAGS) -fdefault-real-8 -fdefault-double-8 
FFLAGS          := $(FFLAGS) -ffree-line-length-none
FFLAGS          := $(FFLAGS) -x f95-cpp-input
## TODO find equivalents of this ifort options
#FFLAGS          := $(FFLAGS)  -xT -align 

FFLAGS_FIXEDFORM = $(FFLAGS) -ffixed-form  -fbounds-check $(FFLAGS_OMP) 
FFLAGS_NOBOUNDS  = $(FFLAGS) -ffixed-form                 $(FFLAGS_OMP) 
FFLAGS          := $(FFLAGS)               $(FFLAGS_OMP) 

#FFLAGS_NO_OMP   = -fpp -r8 -i4 -O2 -xT -align
FFLAGS          := $(FFLAGS) #-DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0
CFLAGS          := -DSTAN_FLAG 
# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 0
USE_MUMPS = 0
USE_WSMP   = 0

#Scotch library
SCOTCH_HOME  = /builds/jorek_tools/scotch
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME = $(HOME)/intel/mumps/MUMPS_4.7.3
LIB_MUMPS  = $(MUMPS_HOME)/lib/libdmumps.a
LIB_ZMUMPS = $(MUMPS_HOME)/lib/libzmumps.a
INC_MUMPS  = $(MUMPS_HOME)/include/
ORDLIB     = $(MUMPS_HOME)/lib/libpord.a

# PASTIX
PASTIX_HOME      = /builds/jorek_tools/pastix/install
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
LIBFFTW   = -L/opt/cluster/fftw/3.2.2/intel-11.1-069/intel-4.0.0.025 -lfftw3

#LIBLAPACK = $(HOME)/intel/lapack/LAPACK/lapack_LINUX.a -lblas
#LIBLAPACK = -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core 
#SCALAP    = -L$(HOME)/intel/scalapack/SCALAPACK -lscalapack
#SCALAP    = -lmkl_scalapack_lp64
#BLACS_HOME = $(HOME)/intel/blacs/BLACS/LIB
#BLACS      = $(BLACS_HOME)/blacs_MPI-LINUX-0.a		\
#	     $(BLACS_HOME)/blacsF77init_MPI-LINUX-0.a 	\
#	     $(BLACS_HOME)/blacsCinit_MPI-LINUX-0.a   	\
#	     $(BLACS_HOME)/blacs_MPI-LINUX-0.a 
#BLACS     = -lmkl_blacs_intelmpi20_lp64

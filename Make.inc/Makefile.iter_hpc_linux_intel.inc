# Configuration file for jorek

# model directory
MODEL = model199

# Fortran compiler 
FC               = mpiifort -mt_mpi
CC               = mpiicc   -mt_mpi
CXX              = mpiicpc  -mt_mpi
FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -msse4.2 -align
FFLAGS          := $(FFLAGS) -O2 -mcmodel=medium -DUSE_BLOCK
#FFLAGS          := $(FFLAGS) -O2 -mcmodel=medium -DUSE_BLOCK -Dfullmhd 


# --- DEBUGGING
#   --- Debug symbols for debuggers
#DEBUGFLAGS  = -g
#   --- Check array boundaries
#DEBUGFLAGS += -check bounds -traceback -fpe-all=0 -ftrapuv
#   --- Generate and check interfaces (__genmod files)
#DEBUGFLAGS += -gen-interfaces -warn interfaces
#   --- Check many other things
#DEBUGFLAGS += -warn all,nounused -check all,noarg_temp_created -debug all -debug-parameters -fstack-security-check -ftrapuv
#   --- Even more checks and warnings
#DEBUGFLAGS += -warn all -check all

FFLAGS_FIXEDFORM := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)
FFLAGS_NOBOUNDS  := $(FFLAGS)               $(FFLAGS_OMP)
FFLAGS           := $(FFLAGS) $(DEBUGFLAGS) $(FFLAGS_OMP)

FFLAGS_NO_OMP     = -fpp -r8 -i4 -O2 -msse4.2 -align 

# Solvers dependencies
USE_HIPS         = 0
USE_PASTIX       = 1
USE_BLOCK        = 1
USE_PASTIX_MURGE = 0
USE_MUMPS        = 1
USE_FFTW         = 1

# Directory containing the user-compiled libraries (e.g., MUMPS, PASTIX, ...)
LIBDIR = $(HOME)/lib

#Scotch library
SCOTCH_HOME ?=  $(LIBDIR)/scotch/scotch_5.1
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lptscotch -lptscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# MUMPS
MUMPS_HOME ?=  $(LIBDIR)/MUMPS_4.8.3
LIB_MUMPS   = -L$(MUMPS_HOME)/lib/ -ldmumps -lmumps_common
LIB_ZMUMPS  = -lzmumps
INC_MUMPS   =  $(MUMPS_HOME)/include/
ORDLIB      = -L$(METIS_HOME)/lib -lpord -lmetis

# PASTIX
PASTIX_HOME     ?=  $(LIBDIR)/pastix/install/
LIB_PASTIX       = -L$(PASTIX_HOME) -lpastix_murge -lpastix $(LIB_SCOTCH)
LIB_PASTIX_MURGE = -L$(PASTIX_HOME) -lpastix_murge -lpastix $(LIB_SCOTCH)
LIB_PASTIX_BLAS  = 
INC_PASTIX       = -I$(PASTIX_HOME) $(INC_SCOTCH)

#HIPS
LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#FFTW
LIBFFTW   = -L$(FFTWDIR)/lib -lfftw3
INC_FFTW  = -I$(FFTWDIR)/include

#LAPACK
MKL       = $(MKLROOT)/lib/intel64
LIBLAPACK = -L$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -limf
SCALAP    = -L$(MKL) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 

#Dierckx
LIBDIERCKX = $(LIBDIR)/Dierckx/dierckx.a

#PPPLIB
PPPLIB = -L$(PPPLIB_HOME)/lib -lppplib

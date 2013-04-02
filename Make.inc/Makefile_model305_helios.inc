# Configuration file for jorek on helios using global libraries
# including flags for model305, model710

# model directory
#MODEL = model305
MODEL = model710

# for model305 with global outputs and alternative current source
#DEFINES = -DJECCD -Daltcs -DUSE_BLOCK

# for full mhd model710
DEFINES = -Dfullmhd -DMULTI_FUN -DUSE_BLOCK

# Fortran compiler 
FC     = mpiifort
FFLAGS_OMP       = -openmp -mt_mpi 
FFLAGS           = -r8 -i4
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -xAVX -align
FFLAGS          := $(FFLAGS) -O2 -mcmodel=medium

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

FFLAGS_NO_OMP    = -fpp -r8 -i4 -O2 -xAVX -align 

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 1
USE_MUMPS = 0

# Directory containing the user-compiled libraries (e.g., MUMPS, PASTIX, ...)
#LIBDIR = /csc/home1/huijsmg/lib

#Scotch library
#SCOTCH_HOME  = /csc/home1/huijsmg/lib/scotch/scotch_5.1
#LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lptscotch -lptscotcherr
#INC_SCOTCH   = -I$(SCOTCH_HOME)/include

SCOTCH_HOME ?= $(SCOTCH_DIR)
LIBSCOTCH   = -L $SCOTCH_DIR/lib -lscotch -lscotcherr $SCOTCH_DIR/lib/libesmumps.a
LIB_SCOTCH   = -L$SCOTCH_DIR/lib -lscotch -lscotcherr $SCOTCH_DIR/lib/libesmumps.a
INC_SCOTCH   = -I$SCOTCH_DIR/include

PASTIX_HOME     ?= $(PASTIX_DIR)
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs`
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`

# PASTIX
#PASTIX_HOME      = /csc/home1/huijsmg/lib/pastix/install
#LIB_PASTIX       = -L$(PASTIX_HOME) -lpastix_murge -lpastix $(LIB_SCOTCH)
#LIB_PASTIX_MURGE = -L$(PASTIX_HOME) -lpastix_murge -lpastix $(LIB_SCOTCH)
#LIB_PASTIX_BLAS  = 
#INC_PASTIX       = -I$(PASTIX_HOME) $(INC_SCOTCH)


#HIPS
LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#FFTW
#LIBFFTW   =-L$(LIBHOME)/fftw-3.2.2/install/lib -lfftw3 

# LIBFFTW
FFTW_HOME     ?= /csc/softs/fftw/fftw-3.3/intel-12.0.5.220/intelmpi-4.0.3/default/
LIBFFTW = $(FFTW_HOME)/lib/libfftw3.a

#LAPACK
#MKL       = /shared/hpc/compiler/intel/composerxe-2011.2.137/mkl/lib/intel64
#LIBLAPACK = -L$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5
#SCALAP    = -L$(MKL) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 
LIBLAPACK = -L$MKLROOT/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5
SCALAP    = -L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64
BLACS     = -L$MKLROOT/lib/intel64 -lmkl_blacs_intelmpi_lp64

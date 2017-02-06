# Configuration file for jorek
# Version for ITER regression testing platform using gfortran compilers

# model directory
MODEL = model199

# Fortran compiler 
FC               = h5pfc
FFLAGS          := $(FFLAGS) -mcmodel=medium -DUSE_BLOCK
CC = h5pcc

# Solvers dependencies
USE_HIPS         = 0
USE_PASTIX       = 1
USE_PASTIX_MURGE = 0
USE_MUMPS        = 1
USE_FFTW         = 1
USE_HDF5         = 1

# Directory containing the user-compiled libraries (e.g., MUMPS, PASTIX, ...)
LIBDIR = $(HOME)/lib

#Scotch library
LIB_SCOTCH   = -lptscotch -lptscotcherr

# MUMPS
LIB_MUMPS   = -ldmumps -lmumps_common
LIB_ZMUMPS  = -lzmumps
ORDLIB      = -lpord -lmetis
INC_MUMPS   = /work/imas/opt/EasyBuild/software/MUMPS/4.10.0-goolf-1.5.16-metis/include # sorry sorry sorry but
# /work/imas/opt/EasyBuild/modules/all/MUMPS/4.10.0-goolf-1.5.16-metis does not define MUMPS_DIR and 
# gfortran does not follow CPATH apparently

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

#LAPACK
MKL       = $(MKLROOT)/lib/intel64
LIBLAPACK = -L$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -limf
SCALAP    = -L$(MKL) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 

#Dierckx
LIBDIERCKX = $(LIBDIR)/Dierckx/dierckx.a

HDF5INCLUDE ?=. # dummy

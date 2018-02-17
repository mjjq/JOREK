# Configuration file for jorek
# Version for ITER regression testing platform using intel compilers
# Not specifically intended for ifort 2017, just a file made in 2017.

# model directory
MODEL = model199

# Fortran compiler 
FC               = mpiifort
FFLAGS          := $(FFLAGS) -mcmodel=medium -DUSE_BLOCK
CC = mpiicc
CXX = mpiicpc

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

#HDF5
HDF5_HOME=/work/imas/opt/hdf5/1.8.15
HDF5INCLUDE=$(HDF5_HOME)/include
HDF5LIB=$(HDF5_HOME)/lib/libhdf5hl_fortran.a \
        $(HDF5_HOME)/lib/libhdf5_hl.a        \
        $(HDF5_HOME)/lib/libhdf5_fortran.a   \
        $(HDF5_HOME)/lib/libhdf5.a           \
	$(ZLIB_HOME)/lib/libz.a

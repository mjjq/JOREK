# Confiduration file for jorek
# model directory
MODEL = model199

# Fortran compiler 
FC     = mpif90 -f90=gfortran
FFLAGS = -g -ffree-line-length-none -x f95-cpp-input -fdefault-real-8 -O0 -fbounds-check -fopenmp # -traceback -xT -align -openmp
FFLAGS_FIXEDFORM = -g -ffixed-form -ffree-line-length-none -x f95-cpp-input -fdefault-real-8 -O0 -fbounds-check -fopenmp # -traceback -xT -align -openmp
FFLAGS_NOBOUNDS = -g -ffixed-form -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -O0  -fopenmp # -traceback -xT -align -openmp
#FFLAGS_NO_OMP = -fpp -r8 -i4 -O2 -xT -align
FFLAGS_NO_OMP = -fopenmp # -traceback -xT -align -openmp
# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 1
USE_MUMPS = 0

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
LIB_PASTIX       = `pastix-conf --libs` 
LIB_PASTIX_MURGE = `pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `pastix-conf --blas`
INC_PASTIX       = `pastix-conf --incs`

#HIPS
LIBHIPS   = -L$(HOME)/hips_trunk/trunk/LIB -lhips -lio $(LIBSCOTCH)
INCHIPS   = -I$(HOME)/hips_trunk/trunk/LIB

#OPENMP
OPENMP_LIB =  -liomp5

# PPPLIB
PPPSRC    = $(HOME)/src/ppplib.f
#intel/jorek/ppplib.a

# LIBFFTW
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

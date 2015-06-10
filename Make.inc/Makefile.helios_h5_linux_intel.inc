# Confiduration file for jorek
# model directory
MODEL = model199

# Fortran compiler 
# with bullxmpi
FC     = mpif90
CC     = mpicc

# With intelmpi
#Matthias# FC     = mpif90 -f90=ifort -mt_mpi
#FC = mpiifort
#CC = mpiicc

# With intelmpi
#FFLAGS_OMP       = -openmp   -mt_mpi
# With bullxmpi
FFLAGS_OMP       = -openmp
FFLAGS           = -r8 -i4
FFLAGS          := $(FFLAGS) -fpp -DUSE_BLOCK
FFLAGS          := $(FFLAGS) -traceback -xAVX -align
# With intelmpi
#FFLAGS          := $(FFLAGS) -O2 #With bullxmpi -DFUNNELED
# With bullxmpi
FFLAGS          := $(FFLAGS) -O2 -DFUNNELED
FFLAGS          := $(FFLAGS) -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0

# --- DEBUGGING
#   --- Debug symbols for debuggers
#DEBUGFLAGS  = -g
#   --- Check array boundaries
#DEBUGFLAGS += -check bounds
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
USE_PASTIX_MURGE = 0
USE_MUMPS = 0
USE_WSMP = 0
USE_HDF5 = 1

#LIBDIR = ~jorek_pstanis/lib
LIBDIR = /csc/home1/jorek_cp/jorek_tools

#Scotch library
SCOTCH_HOME  = $(LIBDIR)/scotch_5.1
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lptscotch -lptscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# PASTIX
# pstanis
#PASTIX_HOME      = $(LIBDIR)/pastix/install
# HOME
PASTIX_HOME      = $(LIBDIR)/pastix_ambiorix/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`

#LAPACK
MKL       = /opt/intel/composer_xe_2013_sp1.4.211/compiler/lib/intel64
LIBLAPACK = -L$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5
SCALAP    = -L$(MKL) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 

#HDF5
HDF5_HOME   = ${HDF5}
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a \
              $(HDF5_HOME)/lib/libhdf5_hl.a \
              $(HDF5_HOME)/lib/libhdf5_fortran.a \
              $(HDF5_HOME)/lib/libhdf5.a -lz



# Confiduration file for jorek
# model directory

# To get ready to compile or execute, please load the following modules :
#       module load intel/14.0.4.211 intelmpi/5.0.1.035

MODEL = model303

# Fortran compiler 
FC     = mpif90 -mt_mpi -openmp
CC     = mpicc  -mt_mpi -openmp
CXX    = mpicxx -mt_mpi -openmp

FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4
FFLAGS          := $(FFLAGS) -fpp -DUSE_BLOCK
FFLAGS          := $(FFLAGS) -traceback -xAVX -align
FFLAGS          := $(FFLAGS) -O2 -DFUNNELED  # Optimization level
#FFLAGS          := $(FFLAGS) -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=0

CFLAGS=-DSTAN_FLAG
# For PaStiX after revision 3945, additionally the following is required (Wiki Matthias)
#FFLAGS := $(FFLAGS) -DFUNNELED -DWORLDWAR2
#FFLAGS := $(FFLAGS) -DWORLDWAR2


# --- DEBUGGING
#   --- Debug symbols for debuggers
DEBUGFLAGS  = -g -traceback
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

#FFLAGS_NO_OMP    = -fpp -r8 -i4 -O2 -xSSSE3 -align

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 0
USE_MUMPS = 0
USE_WSMP   = 0
USE_HDF5 = 0

# PASTIX
PASTIX_HOME = /panfs/panasas/cnt0028/rfm2197/SHARED/JOREK/jorek_tools/pastix_release_3184/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`


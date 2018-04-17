# Configuration file for jorek on poincare machine 
# Environment: module load module load intel/12.1.0 mvapich2/1.8.1 mkl/10.3  hwloc/1.6.2_intel
# Pastix release: freestyle (pastix_5.2.2.11.tar.gz)

# model directory
MODEL = model302

# Fortran compiler 
FC                = mpiifort -mt_mpi
CC                = mpiicc   -mt_mpi
CXX               = mpiicpc  -mt_mpi
FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4 -traceback
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -xSSE4.2 -align -DUSE_BLOCK -DTIMER=4 #-DPLENTY_TIMERS  -DMEMTRACE
#FFLAGS          := $(FFLAGS)  
FFLAGS          := ${FFLAGS}  -DWORLDWAR2  -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=1
FFLAGS          := $(FFLAGS) -O2  # Optimization level
FFLAGS          := $(FFLAGS)
FFLAGS          := $(FFLAGS) 
#FFLAGS          := $(FFLAGS) 
#FFLAGS          := $(FFLAGS) -DNORMTRACE
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

#FFLAGS_NO_OMP    = -fpp -r8 -i4 -O2 -xSSSE3 -align

# Solvers dependencies
USE_HIPS   = 0
USE_PASTIX = 1
USE_PASTIX_MURGE = 1
USE_MUMPS = 0
USE_MURGE = 0

# PASTIX
PASTIX_HOME     =  /gpfshome/mds/staff/glatu/jorek_tools/pastix_freestyle/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`


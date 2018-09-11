# Configuration file for jorek

# model directory
MODEL = model199

#YES -DMURGE_PROD_NO_COMM -DMURGE_USE_SEQUENCE -DMURGE_USE_GETLOCALELEMENTLIST  -DMURGE_INTERFACE_MAJOR_VERSION=1 -DMURGE_INTERFACE_MINOR_VERSION=1
#NO -DMURGE_USE_SEQUENCE_HARM -DMURGE_USE_DUPLICATE_ELEMENT

# Fortran compiler 
FC                = mpiifort -mt_mpi
CC                = mpiicc   -mt_mpi
CXX               = mpiicpc  -mt_mpi
FFLAGS_OMP       = -openmp 
FFLAGS           = -r8 -i4 -traceback
FFLAGS          := $(FFLAGS) -fpp
FFLAGS          := $(FFLAGS) -xSSE4.2 -align -DUSE_BLOCK -DTIMER=4 #-DPLENTY_TIMERS  -DMEMTRACE
#FFLAGS          := $(FFLAGS) 
FFLAGS          := ${FFLAGS} 
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
USE_PASTIX_MURGE = 0
USE_MUMPS = 0
USE_HDF5 = 1

#Scotch library
SCOTCH_HOME ?=  /gpfshome/mds/staff/glatu/jorek_tools/scotch_5.1.12
LIB_SCOTCH   = -L$(SCOTCH_HOME)/lib -lscotch -lscotcherr
INC_SCOTCH   = -I$(SCOTCH_HOME)/include

# PASTIX
PASTIX_HOME     =  /gpfshome/mds/staff/glatu/jorek_tools/pastix_release_3184/install
LIB_PASTIX       = `$(PASTIX_HOME)/pastix-conf --libs` 
LIB_PASTIX_MURGE = `$(PASTIX_HOME)/pastix-conf --libs_murge`
LIB_PASTIX_BLAS  = `$(PASTIX_HOME)/pastix-conf --blas`
INC_PASTIX       = `$(PASTIX_HOME)/pastix-conf --incs`

#HDF5
HDF5_HOME   = /gpfshome/mds/staff/glatu/hdf5-1.8.10
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a \
	      $(HDF5_HOME)/lib/libhdf5_hl.a \
	      $(HDF5_HOME)/lib/libhdf5_fortran.a \
	      $(HDF5_HOME)/lib/libhdf5.a -lz \

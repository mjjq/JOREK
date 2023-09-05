MODEL             = model710
FC                = mpiifort
CC                = mpiicc
CXX               = mpiicpc
 
FFLAGS += -O3 -mcmodel=medium
EXTRA_FLAGS += -lstdc++
COMPILER_FAMILY   = intel
# Never use STRUMPACK and PASTIC at the same time 
USE_PASTIX        = 1
USE_PASTIX_MURGE  = 0
USE_PASTIX6       = 0
USE_STRUMPACK     = 0
USE_MUMPS         = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
USE_DIRECT_CONSTRUCTION=0
USE_BICGSTAB=0
LIBS += -liomp5 -pthread -ldl -lm
 
LIBDIR = /raven/u/ihol/lib/intel
 
ifeq (1, $(USE_PASTIX))
 
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKLROOT}/include
 PASTIX_DIR = $(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12
 
 LIB_PASTIX = -L$(PASTIX_DIR) -lpastix -L$(SCOTCH_DIR)/lib -lscotch -lscotcherr
 INC_PASTIX = -I$(PASTIX_DIR)/include
 
 INC_PASTIX += $(MKLINC)
 LIB_PASTIX += $(MKLLIB)
 
endif
 
ifeq (1, $(USE_MUMPS))
# --- MUMPS
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 MUMPS_DIR = $(LIBDIR)/MUMPS_5.4.1
 LIB_MUMPS = -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPS_DIR)/include
 
 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS_EXTRA += $(MKLINC)
 
endif
 
ifeq (1, $(USE_STRUMPACK))
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 GKLIB_HOME = $(LIBDIR)/GKlib/install
 METIS_HOME = $(LIBDIR)/METIS/install
 PARMETIS_HOME = $(LIBDIR)/ParMETIS/install
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK-7.1.3/install
 
 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack
 
 STRUMPACKINC += -I$(METIS_HOME)/include -I$(GKLIB_HOME)/include -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += -L$(METIS_HOME)/lib -lmetis -L$(GKLIB_HOME)/lib -lGKlib -L$(PARMETIS_HOME)/lib -lparmetis
 
 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)
 
 DEFINES += -DNEWSPK
endif
 
# --- FFTW Library
LIBFFTW           = -L$(FFTW_HOME)/lib -lfftw3
INC_FFTW          = -I$(FFTW_HOME)/include
 
# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
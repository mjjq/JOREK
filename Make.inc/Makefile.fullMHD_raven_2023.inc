# JOREK configuration file for compilation on the TOK cluster at IPP Garching
# see also https://www.jorek.eu/wiki/doku.php?id=draco
 
MODEL             = model710
 
FC                = mpiifort
CC                = mpiicx
#mpiicc
CXX               = mpiicpx
#mpiicpc
 
COMPILER_FAMILY   = intel
 
USE_PASTIX        = 1
USE_STRUMPACK     = 0
USE_MUMPS         = 1
USE_PASTIX_MURGE  = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
 
 
LIBDIR = /raven/u/ihol/lib/intel
 
ifeq (1, $(USE_MUMPS))
 MKLLIB+= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKL_HOME}/include
 
 MUMPSDIR = $(LIBDIR)/MUMPS_5.4.1
 LIB_MUMPS = -L$(MUMPSDIR)/lib -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPSDIR)/include
 
 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS += $(MKLINC)
endif
 
 
ifeq (1, $(USE_PASTIX))
 MKLLIB += -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKL_HOME}/include
 
 PASTIX_DIR = $(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12
 
 LIB_PASTIX = -L$(PASTIX_DIR) -lpastix -L$(SCOTCH_DIR)/lib -lscotch -lscotcherrexit
 INC_PASTIX = -I$(PASTIX_DIR)/include
 
 INC_PASTIX += $(MKLINC)
 LIB_PASTIX += $(MKLLIB)
 
endif
 
ifeq (1, $(USE_STRUMPACK))
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 
 METIS_HOME=$(LIBDIR)/metis_git/install
 PARMETIS_HOME=$(LIBDIR)/ParMETIS_git/install
 GKLIB_HOME=$(LIBDIR)/GKlib_git/install
 
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK-7.1.3/install
 
 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack
 
 STRUMPACKINC += -I$(METIS_HOME)/include -I$(GKLIB_HOME)/include
 STRUMPACKLIB += -L$(METIS_HOME)/lib -lmetis -L$(GKLIB_HOME)/lib -lGKlib
 
 STRUMPACKINC += -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += -L$(PARMETIS_HOME)/lib -lparmetis
 
 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)
 
 DEFINES += -DNEWSPK
endif
 
# --- FFTW Library
#LIBFFTW           = -L$(FFTW_HOME)/lib -lfftw3
LIBFFTW           = $(FFTW_HOME)/lib/libfftw3.a
INC_FFTW          = -I$(FFTW_HOME)/include
 
# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
#HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5LIB     = $(HDF5_HOME)/lib/libhdf5hl_fortran.a  $(HDF5_HOME)/lib/libhdf5_hl.a  $(HDF5_HOME)/lib/libhdf5_fortran.a $(HDF5_HOME)/lib/libhdf5.a -L$(HDF5_HOME)/lib  -lz



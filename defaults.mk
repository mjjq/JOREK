DEBUG ?= 0 # Build in release mode by default

# To use the default compiler flags, set $(COMPILER_FAMILY)
# If you do not want to do this, be sure to set OUTPUT_MODULE_COMMAND to the correct
# one for your compiler in Makefile.inc.

# build in directories to not clutter the repository
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)

#TODO put this somewhere
LIBS = $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)


# Do some guessing to get the compiler family if it is unset
# TODO

# Default flags for gfortran
ifeq ($(COMPILER_FAMILY), gnu)
  FLAGS += -O3 -cpp -fopenmp
  FLAGS += -Wall -Wextra
  FLAGS += -Wcharacter-truncation
  FLAGS += -Wimplicit-interface -Wimplicit-procedure
  FLAGS += -Winteger-division
  FLAGS += -Wintrinsics-std
  FLAGS += -Wsurprising
  #FLAGS += -Wrealloc-lhs-all
  FLAGS += -fimplicit-none
  FLAGS += -flto=4 # link-time optimization with 4 jobs
  FLAGS += -fwhole-program
  # options still to be tested
  #FLAGS += -fexternal-blas
  #FLAGS += -ffast-math # better -Ofast
  #more optimization options
  ifeq ($(DEBUG), 1)
    # Debug flags for gfortran, in ascending order of severity
    FLAGS += -g -Og -ggdb -fno-lto
    FLAGS += -fcheck=all
    FLAGS += -ffpe-trap=invalid,zero,overflow -ftrapv \
#	      -finit-real=nan -finit-int=nan -finit-logical=false
#   FLAGS += -Warray-temporaries -Wconversion-extra
  endif

  OUTPUT_MODULE_COMMAND=-J#no space
endif

# Default flags for intel
ifeq ($(COMPILER_FAMILY), intel)
  FLAGS += -openmp
  FLAGS += -implicitnone
  FLAGS += -cpp
  FLAGS += -warn all
  FLAGS += -align
  FLAGS += -ipo -ipo-jobs4 # like -flto for gfortran, see https://software.intel.com/en-us/node/524765
  ifeq ($(DEBUG), 1)
    # Debug flags for ifort, see http://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
    FLAGS += -O0 -g -traceback
    FLAGS += -check all,noarg_temp_created
    FLAGS += -check bounds
    FLAGS += -check uninit
    FLAGS += -ftrapuv
    FLAGS += -debug all -debug-parameters
    FLAGS += -gen-interfaces -warn-interfaces
    FLAGS += -fstack-security-check
    FLAGS += -fpe0
    FLAGS += -assume ieee_fpe_flags # not sure about this one
  else
    FLAGS += -O3
  endif

  OUTPUT_MODULE_COMMAND=-module #space is important
endif

#TODO identify good default flags for XLF
# Correct preprocessor-defines for IBM XLF Compiler
IBM_DEFINES = `echo $(DEFINES) | sed -e 's/^/-WF,/' -e 's/  */,/g'`

INCLUDES2  := $(INCLUDES) $(DEFINES)
ifdef IBMFC
  DEFINES  := $(DEFINES) -DIBM_MACHINE
  INCLUDES := $(INCLUDES) $(IBM_DEFINES)
else
  INCLUDES := $(INCLUDES) $(DEFINES)
endif
# TODO set the option to output module files to a specific directory for XLF


# Save and load modules from $(MODDIR)
FLAGS := $(FLAGS) -I$(MODDIR) $(OUTPUT_MODULE_COMMAND)$(MODDIR)



# Make any object or module file depend on the makefiles
%.o %.mod: Makefile Makefile.inc
%.o: %.c Makefile Makefile.inc
	@echo "*** $< ***"
	$(CC) $(CFLAGS) -c $< -o $@  $(INCLUDES)

# Make rules for specific files
# This is needed because the file stems must match and we do not really want to recreate the directory structure in $(OBJDIR) and $(DEPDIR)
# Also, it will make everything slightly nicer later
# Template for generating object files from source files
# Touch the .mod file again if it exists (because it is not written if there is no change, and this messes with the make rules)
define F90_O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod:: $(1)%.f90
	$$(FC) $$(FFLAGS) $$(FLAGS) $$(DEFINES) $$(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o
	@test -e $(MODDIR)/$$*.mod && touch $(MODDIR)/$$*.mod || true
endef
# Template for generating dependencies from source file
# TODO: alter sfmakedepend to do the regex stuff below
# rule defined with double colon so that we will not look for %.f90.c etc (terminal rule)
define F90_D_TEMPLATE
$(DEPDIR)/%.d:: $(1)%.f90
	@echo "Generating dependencies for $$<"
	@util/sfmakedepend -f $(DEPDIR)/$$(*F).Td $$<; \
	  grep -F "$$(*F).o: " $(DEPDIR)/$$(*F).Td |\
	  sed -e 's/ hdf5[.]mod//g' -e 's/ mpi[.]mod//g' -e 's/ omp_lib[.]mod//g' \
	  -e 's/ iso_c_binding[.]mod//g' -e 's/ iso_fortran_env[.]mod//g' \
	  -e 's/ fruit[.]mod//g' -e 's/ h5lt[.]mod//g' \
	  -e 's/ mgi_module[.]mod//g' \
	  -e 's/^/$(OBJDIR)\//g' -e 's/ \([^ ][^ ]*\)/ $(MODDIR)\/\1/g' \
	  > $(DEPDIR)/$$(*F).d; rm $(DEPDIR)/$$(*F).Td
endef
# Notes above: hdf5, mpi, omp_lib, iso_c_binding, iso_fortran_env, h5lt are system libraries/intrinsics
# mgi_module is removed here because it is not in all models. Add it again explicitly for model5XX
ifeq ($(MODEL_NUMBER), 500)
.obj/jorek2_main.o: .obj/mgi_module.o .mod/mgi_module.mod
endif
ifeq ($(MODEL_NUMBER), 555)
.obj/jorek2_main.o: .obj/mgi_module.o .mod/mgi_module.mod
endif

# This template defines a program $(file_stem)
# which has prerequisites $(OBJDIR)/$(file_stem).o and as determined by the output of obj_deps.sh
# Since make will try to build the .d files first and re-exec every time
# these will be up-to-date
file_stem=$(notdir $(basename $(1)))
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(shell ./util/obj_deps.sh $(OBJDIR)/$(file_stem).o | sort | uniq)
	$$(FC) $$(FFLAGS) $$(FLAGS) -o $(file_stem) $$^ $$(LIBS)
endef



# Use flags
ifeq (model710, $(MODEL))
  DEFINES  := $(DEFINES) -Dfullmhd
endif

ifeq (1, $(USE_FFTW))
  LIBS     := $(LIBS) $(LIBFFTW)
  DEFINES  := $(DEFINES) -DUSE_FFTW
  INCLUDES := $(INCLUDES) $(INC_FFTW)
endif

ifeq (1, $(USE_PASTIX_MURGE))
  LIBS     := $(LIBS) $(LIB_PASTIX_MURGE) $(LIB_PASTIX_BLAS)
  DEFINES  := $(DEFINES) -DUSE_MURGE
  INCLUDES := $(INCLUDES) $(INC_PASTIX)
endif

ifeq (1, $(USE_PASTIX))
  DEFINES  := $(DEFINES) -DUSE_PASTIX
  ifeq (0, $(USE_PASTIX_MURGE))
    LIBS     := $(LIBS) $(LIB_PASTIX) $(LIB_PASTIX_BLAS)
    INCLUDES := $(INCLUDES) $(INC_PASTIX)
  endif
endif

ifeq (1, $(USE_WSMP))
  DEFINES  := $(DEFINES) -DUSE_WSMP
  LIBS     := $(LIBS) $(LIB_WSMP)
endif

ifeq (1, $(USE_HIPS))
  LIBS := $(LIBS) $(LIBHIPS)
  INCLUDES := $(INCLUDES) $(INCHIPS)
  DEFINES  := $(DEFINES) -DUSE_HIPS
endif

ifeq (1, $(USE_MUMPS))
  LIBS := $(LIBS) $(LIB_MUMPS) $(ORDLIB) $(SCALAP) $(BLACS) $(LIBLAPACK) $(LIBBLAS) $(PPPLIB) $(OPENMP_LIB) $(LIBFFTW)
  INCLUDES := $(INCLUDES) -I$(INC_MUMPS)
  DEFINES := $(DEFINES) -DUSE_MUMPS
endif

ifeq (1, $(USE_HDF5))
  LIBS     := $(LIBS) $(HDF5LIB)
  INCLUDES := $(INCLUDES) -I$(HDF5INCLUDE)
  DEFINES  := $(DEFINES) -DUSE_HDF5
else
  $(warning "USE_HDF5=1 is recommended for input/output")
endif

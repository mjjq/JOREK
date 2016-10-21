DEBUG ?= 0 # Build in release mode by default

# To use the default compiler flags, set $(COMPILER_FAMILY)
# If you do not want to do this, be sure to set OUTPUT_MODULE_COMMAND to the correct
# one for your compiler in Makefile.inc.

# build in directories to not clutter the repository
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)

# Do some guessing to get the compiler family if it is unset
ifeq ($(COMPILER_FAMILY),)
  COMPILER_FAMILY := $(shell $(FC) --version | grep -oi 'intel\|gnu' | tr A-Z a-z)
endif

# Default flags for gfortran
ifeq ($(COMPILER_FAMILY), gnu)
  FLAGS += -cpp
  FLAGS += -Wall -Wextra
  FLAGS += -Wno-tabs
  FLAGS += -ffree-line-length-none
  FLAGS += -fdefault-real-8 -fdefault-double-8
  FLAGS += -Wcharacter-truncation
  FLAGS += -Winteger-division
  FLAGS += -Wintrinsics-std
  FLAGS += -Wsurprising
  FLAGS += -Wno-ampersand
  # options still to be tested
  #FLAGS += -fexternal-blas
  #FLAGS += -ffast-math # better -Ofast
  #more optimization options
  ifeq ($(DEBUG), 1)
    # Debug flags for gfortran, in ascending order of severity
    FLAGS += -g -Og -ggdb
    FLAGS += -fimplicit-none
    FLAGS += -Wimplicit-interface -Wimplicit-procedure
    FLAGS += -fcheck=all
    FLAGS += -Wunused-variable
    FLAGS += -ffpe-trap=invalid,zero,overflow -ftrapv \
	      -finit-real=nan -finit-logical=false
  else
    FLAGS += -O3 -cpp
    FLAGS += -flto=4 # link-time optimization with 4 jobs
    CFLAGS += -flto=4 # also for C routines
  endif

  OUTPUT_MODULE_COMMAND=-J#no space
endif

# Default flags for intel
ifeq ($(COMPILER_FAMILY), intel)
  FLAGS += -fpp
  FLAGS += -warn all
  FLAGS += -warn nounused
  FLAGS += -align
  #FLAGS += -ipo -ipo-jobs4 # like -flto for gfortran, see https://software.intel.com/en-us/node/524765
  # Could take a long time on some machines. Test performance increase first
  ifeq ($(DEBUG), 1)
    # Debug flags for ifort, see http://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
    FLAGS += -O0 -g -traceback
    FLAGS += -check all,noarg_temp_created
    FLAGS += -implicitnone
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
ifdef IBMFC
  DEFINES  := $(IBM_DEFINES) -DIBM_MACHINE
endif
# TODO set the option to output module files to a specific directory for XLF


# Save and load modules from $(MODDIR)
FFLAGS := $(FLAGS) $(FFLAGS) -I$(MODDIR) $(OUTPUT_MODULE_COMMAND)$(MODDIR)



# Make any object or module file depend on the makefiles
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@  $(INCLUDES)

# Make rules for specific files
# This is needed because the file stems must match and we do not really want to recreate the directory structure in $(OBJDIR) and $(DEPDIR)
# Also, it will make everything slightly nicer later
# Template for generating object files from source files
# Touch the .mod file again if it exists (because it is not written if there is no change, and this messes with the make rules)
define O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod:: $(1)%.f90
	$$(FC) $$(FFLAGS) $$(EXTRA_FFLAGS) $$(DEFINES) $$(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o
	@test -e $(MODDIR)/$$*.mod && touch $(MODDIR)/$$*.mod || true

$(OBJDIR)/%.o:: $(1)%.f
	$$(FC) $$(FFLAGS) $$(EXTRA_FFLAGS) -fno-implicit-none $$(DEFINES) $$(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o

$(OBJDIR)/%.o:: $(1)%.c
	$$(CC) $$(CFLAGS) $$(DEFINES) $$(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o

$(OBJDIR)/%.o:: $(1)%.cpp
	$$(CXX) $$(CXXFLAGS) $$(DEFINES) $$(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o
endef
# Template for generating dependencies from source file
define F90_D_TEMPLATE
$(DEPDIR)/%.d: $(1)%.f90
	@echo "Generating dependencies for $$<"
	@util/makedepend $$< $(DIRS) > $(DEPDIR)/$$(*F).d
endef
# First call makedepend for use
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
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(shell ./util/obj_deps $(DEPDIR)/$(file_stem).d)
	$$(FC) $$(FFLAGS) $$(EXTRA_FFLAGS) $$(DEFINES) $$(INCLUDES) -o $(file_stem) $$^ $$(LIBS)
endef



LIBS := $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)
DEFINES += -DJOREK_MODEL=$(MODEL_NUMBER) -DUSE_MPI

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


# Do not check to make these files to speed up and clean -d output
Makefile: ;
Makefile.inc: ;
%.mk: ;
%.f90: ;

.mod/version.h:
	@echo "Generate .mod/version.h"
	@`pwd`/util/version.sh 2>/dev/null > $@
	@echo "#define compile_command '$(FC)'" >> $@
	@echo "#define compile_flags '$(FFLAGS) $(EXTRA_FFLAGS)'" >> $@
	@echo "#define compile_includes '$(INCLUDES)'" >> $@
	@echo "#define compile_defines '$(DEFINES)'" >> $@
	@echo "#define compile_libs '$(LIBS)'" >> $@
	-@echo "#define compile_dir '`pwd`'" >> $@
	-@echo "#define compile_time '`date \"+%F %T\"`'" >> $@
	-@echo "#define compile_user '`whoami`'" >> $@
	-@echo "#define compile_machine '`hostname`'" >> $@
	@echo "#define compile_modules '$(LOADEDMODULES)'" >> $@


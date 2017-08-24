DEBUG ?= 0 # Build in release mode by default

# To use the default compiler flags, set $(COMPILER_FAMILY)
# If you do not want to do this, be sure to set OUTPUT_MODULE_COMMAND to the correct
# one for your compiler in Makefile.inc.

# build in directories to not clutter the repository
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)
INCLUDES+=-I$(MODDIR)

# Detect the compiler vendors (sort to remove duplicates)
F_COMPILER_FAMILY :=$(sort $(shell $(FC) --version | grep -oim 1 'intel\|gcc\|gnu' | tr A-Z a-z | sed 's/gcc/gnu/'))
C_COMPILER_FAMILY :=$(sort $(shell $(CC) --version | grep -oim 1 'intel\|gcc\|gnu' | tr A-Z a-z | sed 's/gcc/gnu/'))
CXX_COMPILER_FAMILY :=$(sort $(shell $(CXX) --version | grep -oim 1 'intel\|gcc\|gnu' | tr A-Z a-z | sed 's/gcc/gnu/'))
ifneq ($(F_COMPILER_FAMILY),$(C_COMPILER_FAMILY))
  $(error "Fortran compiler ($(F_COMPILER_FAMILY)) must be same as C compiler ($(C_COMPILER_FAMILY))")
endif
ifneq ($(CXX_COMPILER_FAMILY),$(C_COMPILER_FAMILY))
  $(error "C++ compiler ($(CXX_COMPILER_FAMILY)) must be same as C compiler ($(C_COMPILER_FAMILY))")
endif
ifeq ($(COMPILER_FAMILY),)
  COMPILER_FAMILY = $(F_COMPILER_FAMILY)
endif
ifneq ($(F_COMPILER_FAMILY),$(COMPILER_FAMILY))
  $(warning "*******************************************************************************************")
  $(warning "Set compiler family ($(COMPILER_FAMILY)) is different from $(FC) ($(F_COMPILER_FAMILY))")
  $(warning "*******************************************************************************************")
endif

# Clear some variables (makefile bug?)
F77FLAGS=
F90FLAGS=

# Default flags for GCC
ifeq ($(COMPILER_FAMILY), gnu)
  FLAGS += -cpp -fopenmp
  FLAGS += -Wall -Wextra
  FLAGS += -Wno-unused-variable
  FFLAGS += -Wintrinsics-std
  FFLAGS += -Wcharacter-truncation
  FFLAGS += -Wsurprising -Wno-tabs
  FFLAGS += -ffree-line-length-none
  F77FLAGS += -fdefault-real-8 -fdefault-double-8
  ifeq ($(DEBUG), 1)
    FLAGS  += -g -Og -ggdb -fno-lto
    FLAGS  += -Wunused-variable
    FFLAGS += -fcheck=all
    FFLAGS += -ffpe-trap=invalid,zero,overflow -ftrapv
    FFLAGS += -finit-real=snan -finit-integer=12345678
    FFLAGS += -Wimplicit-interface -Wimplicit-procedure
    FFLAGS += -Wconversion
    F90FLAGS += -fimplicit-none
  endif

  FFLAGS +=-J$(MODDIR)
endif

# Default flags for intel
ifeq ($(COMPILER_FAMILY), intel)
  COMPILER_MAJOR_VERSION=$(shell $(FC) -V 2>&1 | grep -o "Version [0-9]*" | cut -d' ' -f 2)
  FFLAGS += -align
  ifeq ($(shell test $(COMPILER_MAJOR_VERSION) -ge 15; echo $$?),0)
    FLAGS += -qopenmp
  else
    FLAGS += -openmp
  endif
  FFLAGS += -warn all
  FFLAGS += -warn nointerfaces
  FFLAGS += -warn nounused
  FFLAGS += -fpp
  FFLAGS += -r8
  F77FLAGS += -warn nodeclarations
  ifeq ($(DEBUG), 1)
    # Debug flags for ifort, see http://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
    FLAGS += -O0 -g -traceback
    FLAGS += -ftrapuv
    FLAGS += -debug all -debug-parameters
    FLAGS += -fstack-security-check
    FFLAGS += -check all,noarg_temp_created
    FFLAGS += -check bounds
    FFLAGS += -check uninit
    FFLAGS += -init=snan -init=zero -init=arrays
    FFLAGS += -gen-interfaces -warn-interfaces
    F90FLAGS += -implicitnone
  endif

  FFLAGS +=-module $(MODDIR)
endif

#TODO identify good default flags for XLF
# Correct preprocessor-defines for IBM XLF Compiler
IBM_DEFINES = `echo $(DEFINES) | sed -e 's/^/-WF,/' -e 's/  */,/g'`
ifdef IBMFC
  DEFINES  := $(IBM_DEFINES) -DIBM_MACHINE
endif
# TODO set the option to output module files to a specific directory for XLF

# Make rules for specific files
# This is needed because the file stems must match and we do not really want to recreate the directory structure in $(OBJDIR) and $(DEPDIR)
# Also, it will make everything slightly nicer later
# Template for generating object files from source files
# Touch the .mod file again if it exists (because it is not written if there is no change to the interfaces, and this messes with the make rules)
define O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod:: $(1)%.f90
	$$(FC) $$(FLAGS) $$(FFLAGS) $$(F90FLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o
	@test -e $(MODDIR)/$$*.mod && touch $(MODDIR)/$$*.mod || true

$(OBJDIR)/%.o:: $(1)%.f
	$$(FC) $$(FLAGS) $$(FFLAGS) $$(F77FLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o

$(OBJDIR)/%.o:: $(1)%.c
	$$(CC) $$(FLAGS) $$(CFLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o

$(OBJDIR)/%.o:: $(1)%.cpp
	$$(CXX) $$(FLAGS) $$(CFLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o
endef
# Template for generating dependencies from source file
define F90_D_TEMPLATE
$(DEPDIR)/%.d: $(1)%.f90 | $(MODDIR)/version.h
	@echo "Generating dependencies for $$<"
	@cpp -traditional-cpp -dI $$(DEFINES) $$(INCLUDES) $$< | util/makedepend $$< - $(DIRS) > $(DEPDIR)/$$(*F).d
endef

# This template defines a program $(file_stem)
# which has prerequisites $(OBJDIR)/$(file_stem).o and as determined by the output of obj_deps.sh
# Since make will try to build the .d files first and re-exec every time
# these will be up-to-date
file_stem=$(notdir $(basename $(1)))
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(shell ./util/obj_deps $(DEPDIR)/$(file_stem).d)
	$$(FC) $$(FLAGS) $$(EXTRA_FLAGS) $$(DEFINES) $$(INCLUDES) -o $(file_stem) $$^ $$(LIBS)
endef



LIBS += $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)
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

ifeq (1, $(USE_BLOCK))
  DEFINES := $(DEFINES) -DUSE_BLOCK
endif


# Do not check to make these files to speed up and clean -d output
Makefile: ;
Makefile.inc: ;
%.mk: ;
%.f90: ;

.mod/version.h:
	@echo "Generate .mod/version.h"
	@echo "#define RCS_VERSION '`git describe --always --dirty --abbrev 2> /dev/null`'" > $@
	@echo "#define compile_command '$(FC)'" >> $@
	@echo "#define compile_flags '$(FLAGS) $(FFLAGS) $(F90FLAGS) $(EXTRA_FLAGS)'" >> $@
	@echo "#define compile_includes '$(INCLUDES)'" >> $@
	@echo "#define compile_defines '$(DEFINES)'" >> $@
	@echo "#define compile_libs '$(LIBS)'" >> $@
	-@echo "#define compile_dir '`pwd`'" >> $@
	-@echo "#define compile_time '`date \"+%F %T\"`'" >> $@
	-@echo "#define compile_user '`whoami`'" >> $@
	-@echo "#define compile_machine '`hostname`'" >> $@
	@echo "#define compile_modules '$(LOADEDMODULES)'" >> $@


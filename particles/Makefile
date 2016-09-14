# This makefile automatically searches and compiles fortran programs
# Known limitations:
# - All dependencies must be in modules
# - Module name must be equal to file name

# Include jorek-specific things and settings
include Makefile.inc
ifeq (1, $(USE_MUMPS))
  LIBS     := $(LIBS) $(LIB_MUMPS) $(ORDLIB)
  INCLUDES := $(INCLUDES) -I$(INC_MUMPS)
else
  $(warning "USE_MUMPS=1 is required for ADAS support")
endif
ifeq (1, $(USE_HDF5))
  LIBS     := $(LIBS) $(HDF5LIB)
  INCLUDES := $(INCLUDES) -I$(HDF5INCLUDE)
else
  $(warning "USE_HDF5=1 is required for input/output")
endif

# Build in release mode by default
DEBUG ?= 0
FLAGS :=
FLAGS += -DUSE_MPI -Ilibs/flogging/include # For flogging

# Default flags for gfortran
ifeq ($(COMPILER), gnu)
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
  # options still to be tested
  #FFLAGS += -fexternal-blas
  #FFLAGS += -ffast-math # or -Ofast
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
ifeq ($(COMPILER), intel)
  FLAGS += -openmp
  FLAGS += -implicitnone
  FLAGS += -cpp 
  FLAGS += -warn all
  FLAGS += -xAVX -align
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


# build-private directories
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)


# Prepend the default flags so FFLAGS can override it
FFLAGS := $(FLAGS) $(FFLAGS) -I$(MODDIR) $(OUTPUT_MODULE_COMMAND)$(MODDIR)


# default action
all: particles
clean:
	rm -rf $(OBJDIR)
	rm -f tests/fruit_test_driver.f90
cleanall: clean
	rm -rf $(DEPDIR) $(MODDIR)
	find . -name '*.mod' -delete
.PHONY: clean cleanall test fruit docs
test: fruit
fruit:
	PYTHONPATH=$(FRUITPYDIR) python tests/fruit.py
docs: pusher_test compare_mc_coronal
	./pusher_test
	./compare_mc_coronal
	ford particles.md --no-search

# Most of the sources are actually fortran 2003/8 but use the .f90 suffix
sources:=$(shell find libs src tests -iname '*.f90')
depends:=$(patsubst %.f90,$(DEPDIR)/%.d,$(notdir $(sources)))

# Make rules for all directories that we have found sources in
# This is needed because the file stems must match and we do not really want to recreate the directory structure in $(OBJDIR) and $(DEPDIR)
# Also, it will make everything slightly nicer later
# Template for generating object files from source files
define F90_O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod: $(1)%.f90
	$$(FC) $$(FFLAGS) $$(DEFINES) $$(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o
endef
# Template for generating dependencies from source file
# TODO: fork sfmakedepend to do all of this automatically
define F90_D_TEMPLATE
$(DEPDIR)/%.d: $(1)%.f90
	@echo "Generating dependencies for $$<"
	@libs/sfmakedepend/sfmakedepend -f $(DEPDIR)/$$(*F).Td $$<; \
	  grep -F "$$(*F).o: " $(DEPDIR)/$$(*F).Td |\
	  sed -e 's/ hdf5[.]mod//g' -e 's/ mpi[.]mod//g' -e 's/ omp_lib[.]mod//g' \
	  -e 's/ iso_c_binding[.]mod//g' -e 's/ iso_fortran_env[.]mod//g' \
	  -e 's/ fruit[.]mod//g' \
	  -e 's/^/$(OBJDIR)\//g' -e 's/ \([^ ][^ ]*\)/ $(MODDIR)\/\1/g' \
	  > $(DEPDIR)/$$(*F).d; rm $(DEPDIR)/$$(*F).Td
endef
# For each directory containing source files add the templates above
$(foreach dir,$(sort $(dir $(sources))),$(eval $(call F90_O_TEMPLATE,$(dir))))
$(foreach dir,$(sort $(dir $(sources))),$(eval $(call F90_D_TEMPLATE,$(dir))))

# what does this do again?
%.d: ;
-include $(depends)


# This template defines a program $(file_stem)
# which has prerequisites $(OBJDIR)/$(file_stem).o and as determined by the output of obj_deps.sh
# Since make will try to build the .d files first and re-exec every time
# these will be up-to-date
file_stem=$(notdir $(basename $(1)))
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(shell ./util/obj_deps.sh $(OBJDIR)/$(file_stem).o | sort | uniq)
	$$(FC) $$(FFLAGS) -o $(file_stem) $$^ $$(LIBS)
endef

# All programs to compile
PROGRAM_SOURCES = src/particles.f90 \
		  src/static_field_example.f90 \
		  src/jorek_static_fields_boris.f90 \
		  tests/fruit_test_driver.f90 \
		  tests/pushers/pusher_test.f90 \
		  tests/openadas/compare_mc_coronal.f90

$(foreach prog,$(PROGRAM_SOURCES),$(eval $(call PROGRAM_TEMPLATE,$(prog))))

# non-standard object files (extra include options, target-specific variable)
fruit_test_driver: INCLUDES += -I$(INC_FRUIT)
fruit_test_driver: LIBS += -L$(LIB_FRUIT) -lfruit -Wl,-rpath=$(LIB_FRUIT)

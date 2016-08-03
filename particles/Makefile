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
  $(warning "USE_MUMPS=1 is mandatory")
endif
ifeq (1, $(USE_HDF5))
  LIBS     := $(LIBS) $(HDF5LIB)
  INCLUDES := $(INCLUDES) -I$(HDF5INCLUDE)
else
  $(warning "USE_HDF5=1 is mandatory")
endif

# Build in release mode by default
DEBUG ?= 0
FLAGS :=
# Default flags for gfortran
ifeq ($(COMPILER), gnu)
  FLAGS += -O3 -cpp -fopenmp
  FLAGS += -Wall -Wextra
  FLAGS += -Wcharacter-truncation
  FLAGS += -Wimplicit-interface -Wimplicit-procedure
  FLAGS += -Winteger-division
  FLAGS += -Wintrinsics-std
  FLAGS += -Wsurprising
  FLAGS += -Wrealloc-lhs-all
  FLAGS += -fimplicit-none
  FLAGS += -flto=4 # link-time optimization with 4 jobs
  # options still to be tested
  #FFLAGS += -fexternal-blas
  #FFLAGS += -ffast-math
  #more optimization options
  ifeq ($(DEBUG), 1)
    # Debug flags for gfortran, in ascending order of severity
    FLAGS += -g -Og -ggdb -fno-lto
    FLAGS += -fcheck=all
    FLAGS += -ffpe-trap=all -ftrapv -pedantic \
#	      -finit-real=nan -finit-int=nan -finit-logical=false
#   FLAGS += -Warray-temporaries -Wconversion-extra
  endif

  OUTPUT_MODULE_COMMAND=-J
endif

# Default flags for intel
ifeq ($(COMPILER), intel)
  ifeq ($(DEBUG), 1)
    # Debug flags for ifort
    FLAGS += 
  endif

  OUTPUT_MODULE_COMMAND=-module
endif


# build-private directories
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)


# Prepend the default flags so FFLAGS can override it
FFLAGS := $(FLAGS) $(FFLAGS) -I$(MODDIR) $(OUTPUT_MODULE_COMMAND)$(MODDIR)
# OUTPUT_MODULE_COMMAND is -J for gfortran and -module for ifort


# default action
all: test_base_hdf5_io
clean:
	rm -rf $(OBJDIR)
cleanall:
	rm -rf $(DEPDIR) $(MODDIR) $(OBJDIR)
.PHONY: clean cleanall test fruit
test: fruit
fruit:
	PYTHONPATH=$(FRUITPYDIR) python fruit.py



# Most of the sources are actually fortran 2003/8 but use the .f90 suffix
sources:=$(shell find . -iname '*.f90')
depends:=$(patsubst %.f90,$(DEPDIR)/%.d,$(notdir $(sources)))

# Make rules for all directories that we have found sources in
# This is needed because the file stems must match and we do not really want to recreate the directory structure in $(OBJDIR) and $(DEPDIR)
# Also, it will make everything slightly nicer later
# Template for generating object files from source files
define F90_O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod: $(1)%.f90
	$(FC) $(FFLAGS) $(DEFINES) $(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o
endef
# Template for generating dependencies from source file
define F90_D_TEMPLATE
$(DEPDIR)/%.d: $(1)%.f90
	@echo "Generating dependencies for $$<"
	@libs/sfmakedepend/sfmakedepend -f $(DEPDIR)/$$(*F).Td $$<; \
	  grep -F "$$(*F).o: " $(DEPDIR)/$$(*F).Td |\
	  sed -e 's/hdf5[.]mod//g' -e 's/mpi[.]mod//g' -e 's/omp_lib[.]mod//g' \
	  -e 's/iso_c_binding[.]mod//g' -e 's/iso_fortran_env[.]mod//g' \
	  -e 's/fruit[.]mod//g' \
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
# which has prerequisites $(OBJDIR)/$(file_stem).o and as determined by the output of the obj_deps.sh script
# Since make will try to build the .d files first and re-exec every time
# these will be up-to-date
file_stem=$(notdir $(basename $(1)))
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(shell ./util/obj_deps.sh $(OBJDIR)/$(file_stem).o | sort | uniq)
	$(value FC) $(value FFLAGS) $(2) -o $(file_stem) $$^ $(value LIBS)
endef

# All programs to compile
PROGRAM_SOURCES = tests/input_output/test_base_hdf5_io

$(foreach prog,$(PROGRAM_SOURCES),$(eval $(call PROGRAM_TEMPLATE,$(prog))))

# non-standard programs (extra linking options)
FRUIT_LINK_OPTIONS=-L$(LIB_FRUIT) -lfruit -Wl,-rpath=$(LIB_FRUIT)
$(eval $(call PROGRAM_TEMPLATE,fruit_test_driver,$(FRUIT_LINK_OPTIONS)))

# non-standard object files (extra include options)
$(OBJDIR)/fruit_test_driver.o $(MODDIR)/fruit_test_driver.mod: tests/fruit_test_driver.f90
	$(FC) $(FFLAGS) $(DEFINES) $(INCLUDES) -I$(INC_FRUIT) -c $< -o $*.o

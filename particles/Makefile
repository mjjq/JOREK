# Include jorek-specific things and settings
include Makefile.inc
ifeq (1, $(USE_MUMPS))
  LIBS := $(LIBS) $(LIB_MUMPS) $(ORDLIB) $(SCALAP) $(BLACS) $(LIBLAPACK) $(LIBBLAS) $(PPPLIB) $(OPENMP_LIB) $(LIBFFTW)
  INCLUDES := $(INCLUDES) -I$(INC_MUMPS)
  DEFINES := $(DEFINES) -DUSE_MUMPS
endif
ifeq (1, $(USE_HDF5))
  LIBS     := $(LIBS) $(HDF5LIB)
  INCLUDES := $(INCLUDES) -I$(HDF5INCLUDE)
  DEFINES  := $(DEFINES) -DUSE_HDF5
endif

# TODO clean up FFLAGS from JOREK
#FFLAGS += -std=f2008ts #TODO set standard and remove old style code

# output directories
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)

FFLAGS += -I$(MODDIR) $(OUTPUT_MODULE_COMMAND)$(MODDIR)
# OUTPUT_MODULE_COMMAND is -J for gfortran and -module for ifort

# path of outfile of a specific source file
df = $(DEPDIR)/$(*F)
of = $(OBJDIR)/$(*F)
mf = $(MODDIR)/$(*F)


# Make dependency information
#MAKEDEPEND = gfortran -cpp -MM -MG -MT '$(patsubst %.f90,%.o,$<)' $< > $(df).Td
# This transform takes the files outputted by sfmakedepend
# and for every .mod it adds $(MODDIR)/
define MAKEDEPEND_TRANSFORM
cp $(df).Td $(df).d; \
sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
-e '/^$$/ d' -e 's/$$/ :/' < $(df).Td >> $(df).d; \
rm -f $(df).Td
endef


all: test_base_hdf5_io
clean:
	rm -rf $(OBJDIR)
cleanall:
	rm -rf $(DEPDIR) $(MODDIR) $(OBJDIR)

.PHONY: clean cleanall



# Most of the sources are actually fortran 2003 but use the .f90 suffix
# TODO: rename all source files to .f03 and update Makefile
sources:=$(shell find . -iname '*.f90')
depends:=$(patsubst %.f90,$(DEPDIR)/%.d,$(notdir $(sources)))

# Make rules for all directories that we have found sources in
define F90_O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod: $(1)%.f90 $(DEPDIR)/%.d
	$(FC) $(FFLAGS) $(DEFINES) $(INCLUDES) -c $$< -o $(OBJDIR)/$$*.o
endef
define F90_D_TEMPLATE
$(DEPDIR)/%.d: $(1)%.f90
	@echo "Generating dependencies for $$<"
	@libs/sfmakedepend/sfmakedepend -f $(DEPDIR)/$$(*F).d $$<; \
	  sed -i $(DEPDIR)/$$(*F).d -e '/^#.*/d' \
	  -e '/.mod:/d' \
	  -e 's/hdf5.mod//g' \
	  -e 's/mpi.mod//g' \
	  -e 's/omp_lib.mod//g' \
	  -e 's/iso_c_binding.mod//g' \
	  -e 's/iso_fortran_env.mod//g' \
	  -e 's/^/$(OBJDIR)\//g' \
	  -e 's/ \([^ ][^ ]*\)/ $(MODDIR)\/\1/g'
endef
$(foreach dir,$(sort $(dir $(sources))),$(eval $(call F90_O_TEMPLATE,$(dir))))
$(foreach dir,$(sort $(dir $(sources))),$(eval $(call F90_D_TEMPLATE,$(dir))))


%.d: ;
-include $(depends)




object_dependencies=`./util/obj_deps.sh $(OBJDIR)/$(notdir $(basename $(1))).o`

#$(info $(call all_object_dependencies,test_base_hdf5_io))
file_stem=$(notdir $(basename $(1)))
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(DEPDIR)/$(file_stem).d
	$(value FC) $(value FFLAGS) $(value INCLUDES) -c $(1) -o $(OBJDIR)/$(file_stem).o
	$(value FC) $(value FFLAGS) -o $(file_stem) $(value LIBS) $(object_dependencies)
endef


PROGRAM_SOURCES = tests/input_output/test_base_hdf5_io.f90

# Create all standard make rules
$(foreach prog,$(PROGRAM_SOURCES),$(eval $(call PROGRAM_TEMPLATE,$(prog))))

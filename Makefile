include Makefile.inc

SED   ?= sed # gnu-sed command
AWK   ?= awk # gnu-awk command
IBMFC ?=     # IBM compiler flag (for FORTRAN symbol defs)

JOREK_DIR = `pwd`

DIRS =  timing				\
        models/$(MODEL) 		\
        datatypes			\
	models				\
	communication			\
	elements			\
	grids				\
	grids/grid_utils		\
	matrix				\
	solvers 			\
	plots				\
	diagnostics			\
	diagnostics/new_diag		\
	particles                       \
	refinement			\
	postproc			\
	tools                           \
	vacuum				


LIBS = $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)

NODEPS = clean cleanall cleandep forcheck forcheck_poincare \
    forcheck_rst_bin2hdf5 forcheck_rst_hdf52bin\
    forcheck_four forcheck_postproc forcheck_connection2 forcheck_jorek2vtk	\
    forcheck_fieldlines_vtk forcheck_jorek2vtk_3d forcheck_diagno		\
    forcheck_strikes forcheck_jorek_to_helena forcheck_import_eqdsk            

# If we have neither HIPS or PASTIX_MURGE we need to
# Use fake murge.

INCLUDES  =   -I. $(patsubst %,-I%/,$(DIRS)) $(INCMURGE)

VPATH = $(MAIN_MODEL_DIR) $(MODEL_DIR) $(DATATYPES_DIR) $(SOLVERS_DIR)

.SUFFIXES: .o .f90 .f .c

MODEL_NUMBER = `echo $(MODEL) | sed -e 's/model//'`
DEFINES := $(DEFINES) -DJOREK_MODEL=$(MODEL_NUMBER)

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
endif

# Correct preprocessor-defines for IBM XLF Compiler
IBM_DEFINES = `echo $(DEFINES) | sed -e 's/^/-WF,/' -e 's/  */,/g'`

INCLUDES2  := $(INCLUDES) $(DEFINES)
ifdef IBMFC
  DEFINES  := $(DEFINES) -DIBM_MACHINE
  INCLUDES := $(INCLUDES) $(IBM_DEFINES)
else
  INCLUDES := $(INCLUDES) $(DEFINES)
endif

# include the description for each module
include $(patsubst %,%/module.mk,$(DIRS))

PROGRAMS := JOREK2_MAIN JOREK2_POINCARE RST_BIN2HDF5 RST_HDF52BIN           \
            JOREK2_POSTPROC JOREK2_POVRAY JOREK2_CONNECTION2 JOREK2_STRIKES \
            ENBIGGEN JORDEL JORPOL JOREK2VTK JOREK2FLVTK JOREK2VTK3D        \
            JOREK2_FOUR JOREK_EXTRACT_DATA JOREK2_DIAGNO JOREK_TO_HELENA    \
            JOREK2_TARGET2VTK JOREK2_POWERS JOREK2_IMPORT_PERTURBATION      \
            NEW_DIAG_DEMO PENNING_TEST JOREK2_PARTICLES SIMON_PARTICLE_TEST \
	    PARTICLE_TEST PROJECT_PARTICLES_VTK COUNT_PARTICLES_VTK 	    \
	    DUMP_PARTICLES_VTK PARTICLE_FLUX_COORDINATES PARTICLE_FLUX_COORDINATE_DIFFUSION

# Add the common sources to all these programs
$(foreach prog,$(PROGRAMS),$(eval $(prog)_SRC += $(ALL_BINARIES_SRC) $(PPPSRC)))

# Add the jorek2_main sources to the simon particle test (needs to calculate equilibrium)
SIMON_PARTICLE_TEST_SRC += $(JOREK2_MAIN_SRC)
SIMON_PARTICLE_TEST_SRC := $(sort $(SIMON_PARTICLE_TEST_SRC))

# Add extra source files
JOREK2_MAIN_SRC += jorek2_main.f90

# sort source dependencies
SRC_DEP := $(sort $(foreach prog,$(PROGRAMS),$(call $(prog)_SRC)))
SRC_DEP := $(filter %.f90, $(SRC_DEP)) $(filter %.f, $(SRC_DEP)) diagnostics/hdf5_library.important

# Create $(prog)_OBJ files by replacing $(SOURCE_SUFFIXES) -> .o
SOURCE_SUFFIXES = .f90 .f .c
all_src_to_obj = $(foreach SUFFIX,$(SOURCE_SUFFIXES), $(patsubst %$(SUFFIX),%.o,$(filter %$(SUFFIX),$(1))))
$(foreach prog,$(PROGRAMS),$(eval $(prog)_OBJ := $(call all_src_to_obj,$(call $(prog)_SRC))))
# equivalent to JOREK2_MAIN_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_MAIN_SRC))) \ (and for .f and .c)


MOD_FILES=`find . -name "*.mod"`
MAIN = jorek_$(MODEL)

.PHONY: version

all: $(MAIN) version

cleanall : clean cleandep cleangenmod
	@echo ">> Deleting some executables"
	-@rm -f rst_bin2hdf5 rst_hdf52bin

clean :	
	@echo ">> Deleting Object Files <<"
	-@rm -f $(JOREK2_MAIN_OBJ) $(JOREK2_FOUR_OBJ) $(JOREK2_POSTPROC_OBJ)
	@echo ">> Deleting Module Files <<"
	-@rm -f $(MOD_FILES)

cleandep:
	@echo ">> Deleting Dependency Files <<"
	-@rm -f *.dep */*.dep */*/*.dep

cleangenmod:
	@echo ">> Deleting auto-generated interface files (*__genmod*) if any <<"
	-@rm -f *__genmod*

version.h: version

version:
	@echo "Generate version.h"
	@rm -f version.h
	@$(JOREK_DIR)/util/version.sh  2>/dev/null
	@echo "#define compile_command '$(FC)'" >> version.h
	@echo "#define compile_flags '$(FFLAGS)'" >> version.h
	@echo "#define compile_includes '$(INCLUDES)'" >> version.h
	@echo "#define compile_defines '$(DEFINES)'" >> version.h
	@echo "#define compile_libs '$(LIBS)'" >> version.h
	-@echo "#define compile_dir '`pwd`'" >> version.h
	-@echo "#define compile_time '`date \"+%F %T\"`'" >> version.h
	-@echo "#define compile_user '`whoami`'" >> version.h
	-@echo "#define compile_machine '`hostname`'" >> version.h
	@echo "#define compile_modules '$(LOADEDMODULES)'" >> version.h

%.dep:%.f90
	@echo "Generating Dependencies for$(patsubst %.f90, %.o, $<)"
	@cpp $(INCLUDES2) < $< 2>/dev/null| grep -i "^[[:space:]]*use " | grep -v -i "iso_c_binding" | $(SED) 's/\,.*//' | 					\
		$(AWK) -v file_o=" $(patsubst %.f90, %.o, $<)" '{print file_o" : "tolower($$2)".mod"}' >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES2) < $< 2>/dev/null| grep -i "^[[:space:]]*include " | grep -i "f90" |                                                   \
		$(SED) "s|^[ \t]*include[ \t]*['\"]\([^!'\" \t]*.f90\)['\"].*|`echo $< | $(SED) 's|.f90||'`.o : `echo $< | $(SED) 's|\(.*/\).*|\1|'`\1|I"        \
		>> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES2) < $< 2>/dev/null|											\
		$(SED) -n 's@^\# *[0-9][0-9]* *"\([^"]*\)".*@$(patsubst %.f90, %.o, $<): \1@p' | 				\
		sort | uniq | grep -v ": /" | grep -v ": <">> $@.tmp || touch $@.tmp;
	@grep -q -i "^[[:space:]]*module" $< ; 											\
	if [ $$? -eq 0 ]; then 													\
		grep -i "^[[:space:]]*module" $<										\
		| sed 's/\r$$//' | $(AWK) -v file="$(patsubst %.f90, %.o, $<)" '{print tolower($$2)".mod : "file}' >> $@.tmp || touch $@.tmp;	\
	fi;
	-@$(SED) -e "s/murge.inc//g" -e "s/dmumps_struc.h//g" < $@.tmp > $@ || touch $@
	-@rm -f $@.tmp

%.dep: %.f
	@echo "Generating Dependencies for$(patsubst %.f, %.o, $<)"
	@cpp $(INCLUDES) < $< 2>/dev/null| grep -i "^[[:space:]]*use " | 							\
		$(AWK) -v file_o="$(patsubst %.f, %.o, $<)" '{print file_o" : "tolower($$2)".mod"}' >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES) < $< 2>/dev/null| 											\
		$(SED) -n "s@[ ]+include[ ]+'\([^']*\)*.*@$(patsubst %.f, %.o, $<) : \1@pi" >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES) < $< 2>/dev/null| 											\
		$(SED) -n 's@[ ]+include[ ]+"\([^"]*\)*.*@$(patsubst %.f, %.o, $<) : \1@pi' >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES) < $< 2>/dev/null|											\
		$(SED) -n 's@^\# *[0-9][0-9]* *"\([^"]*\)".*@$(patsubst %.f, %.o, $<): \1@p' | 					\
		sort | uniq | grep -v ": /" | grep -v ": <">> $@.tmp || touch $@.tmp;
	@grep -q -i "^[[:space:]]*module" $< ; 											\
	if [ $$? -eq 0 ]; then 													\
		grep -i "^[[:space:]]*module" $<										\
		| sed 's/\r$$//' | $(AWK) -v file="$(patsubst %.f, %.o, $<)" '{print tolower($$2)".mod : "file}' >> $@.tmp || touch $@.tmp;	\
	fi;
	-@$(SED) -e "s/murge.inc//g" -e "s/dmumps_struc.h//g" < $@.tmp > $@ || touch $@
	-@rm -f $@.tmp


# This line defines the default program template (after https://www.gnu.org/software/make/manual/html_node/Eval-Function.html)
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): version $(1) $$($(shell echo $(notdir $(basename $(1))) | tr '[:lower:]' '[:upper:]')_OBJ)
	$(value FC) $(value FFLAGS) $(value INCLUDES) -c $(1) -o $(call all_src_to_obj,$(1))
	$(value FC) $(value FFLAGS_OMP) $(call all_src_to_obj,$(1)) $$($(shell echo $(notdir $(basename $(1))) | tr '[:lower:]' '[:upper:]')_OBJ) \
	-o $(JOREK_DIR)/$(notdir $(basename $(1))) $(value LIBS)
endef

# Below is an example:
#jorek2_poincare : version diagnostics/jorek2_poincare.f90 $(JOREK2_POINCARE_OBJ) 
#	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_poincare.f90 -o diagnostics/jorek2_poincare.o
#	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_poincare.o $(JOREK2_POINCARE_OBJ) \
#	 -o $(JOREK_DIR)/jorek2_poincare  $(LIBS)

PROGRAM_SOURCES = diagnostics/jorek2_poincare.f90       \
		  diagnostics/rst_bin2hdf5.f90          \
		  diagnostics/rst_hdf52bin.f90          \
		  diagnostics/jorek2_four.f90           \
		  postproc/jorek2_postproc.f90          \
		  diagnostics/jorek2_povray.f90         \
		  diagnostics/jorek2_connection.f90     \
		  diagnostics/enbiggen.f90              \
		  diagnostics/jordel.f90                \
		  diagnostics/jorpol.f90                \
		  diagnostics/jorek2vtk.f90             \
		  diagnostics/jorek2_diagno.f90         \
		  diagnostics/jorek_to_helena.f90       \
		  diagnostics/jorek2_target2vtk.f90     \
		  diagnostics/jorek2_powers.f90         \
		  diagnostics/new_diag_demo.f90 	\
		  particles/jorek2_particles.f90 	\
		  particles/count_particles_vtk.f90 	\
		  particles/dump_particles_vtk.f90 	\
		  particles/particle_flux_coordinates.f90 \
		  particles/particle_flux_coordinate_diffusion.f90 \
		  particles/project_particles_vtk.f90

# Create all standard make rules
$(foreach prog,$(PROGRAM_SOURCES),$(eval $(call PROGRAM_TEMPLATE,$(prog))))


# The below are non-standard make rules
$(MAIN) : version $(JOREK2_MAIN_OBJ)
	$(FC) $(FFLAGS_OMP) \
	$(JOREK2_MAIN_OBJ) \
	 -o $(MAIN) $(INCLUDES) $(LIBS)

jorek_extract_data : version $(JOREK_EXTRACT_DATA_OBJ) diagnostics/jorek_extract_data.o
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek_extract_data.o \
	$(JOREK_EXTRACT_DATA_OBJ)		\
	 -o $(JOREK_DIR)/jorek_extract_data $(INCLUDES) $(LIBS) $(LIBFFTW)

jorek2_connection_stan : version diagnostics/jorek2_connection_stan.f90 $(JOREK2_CONNECTION2_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_connection_stan.f90 -o diagnostics/jorek2_connection_stan.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_connection_stan.o $(JOREK2_CONNECTION2_OBJ) \
	 -o $(JOREK_DIR)/jorek2_connection_stan $(LIBS)

jorek2_strikes : version diagnostics/jorek2_strikes_ordered.f90 $(JOREK2_STRIKES_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_strikes_ordered.f90 -o diagnostics/jorek2_strikes_ordered.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_strikes_ordered.o $(JOREK2_STRIKES_OBJ) \
	 -o $(JOREK_DIR)/jorek2_strikes $(LIBS)

jorek2_fieldlines_vtk : version diagnostics/jorek2_fieldlines_vtk.f90 $(JOREK2FLVTK_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_fieldlines_vtk.f90 -o diagnostics/jorek2_fieldlines_vtk.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_fieldlines_vtk.o $(JOREK2FLVTK_OBJ) \
	 -o $(JOREK_DIR)/jorek2_fieldlines_vtk $(LIBS)

jorek2vtk_3d : version diagnostics/jorek2vtk_3d.f90 $(JOREK2VTK3D_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2vtk_3d.f90 -o diagnostics/jorek2vtk_3d.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2vtk_3d.o $(JOREK2VTK3D_OBJ) \
	 -o $(JOREK_DIR)/jorek2vtk_3d $(LIBS)

jorek2vtk_GaussVortTerms : diagnostics/jorek2vtk_GaussVortTerms.f90 $(JOREK2VTK_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2vtk_GaussVortTerms.f90 -o diagnostics/jorek2vtk_GaussVortTerms.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2vtk_GaussVortTerms.o $(JOREK2VTK_OBJ) \
	 -o $(JOREK_DIR)/jorek2vtk_GaussVortTerms $(LIBS)

jorek2_stan : version diagnostics/jorek2_stan.f90 $(JOREK2VTK_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_stan.f90 -o diagnostics/jorek2_stan.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_stan.o $(JOREK2VTK_OBJ) \
	 -o $(JOREK_DIR)/jorek2_stan $(LIBS)

import_eqdsk : version util/import_eqdsk.f90
	$(FC) -c util/import_eqdsk.f90 -o util/import_eqdsk.o
	$(FC)  util/import_eqdsk.o $(JOREK_DIR)/import_eqdsk $(LIBS)

jorek2_import_perturbation : version jorek2_import_perturbation_new_flags
jorek2_import_perturbation_new_flags: version FFLAGS += -DIMPORT_PERTURBATIONS
jorek2_import_perturbation_new_flags: version jorek2_import_perturbation_tmp

jorek2_import_perturbation_tmp : version diagnostics/jorek2_import_perturbation.f90 $(JOREK2_IMPORT_PERTURBATION_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_import_perturbation.f90 -o diagnostics/jorek2_import_perturbation.o
	$(FC) $(FFLAGS) diagnostics/jorek2_import_perturbation.o $(JOREK2_IMPORT_PERTURBATION_OBJ) \
	-o $(JOREK_DIR)/jorek2_import_perturbation $(LIBS)

penning_test : non_regression_tests/standalone/penning_test/penning.f90 $(PENNING_TEST_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) non_regression_tests/standalone/penning_test/penning.f90 \
	$(PENNING_TEST_OBJ) -o $(JOREK_DIR)/penning_test \
	$(INCLUDES) $(LIBS)

simon_particle_test : non_regression_tests/standalone/simon_particle_test/simon_particle_test.f90 $(SIMON_PARTICLE_TEST_OBJ)
	$(FC) $(FFLAGS) non_regression_tests/standalone/simon_particle_test/simon_particle_test.f90 \
	$(SIMON_PARTICLE_TEST_OBJ) -o $(JOREK_DIR)/simon_particle_test \
	$(INCLUDES) $(LIBS)

include all_rules.mk
ifeq (0, $(words $(foreach word, ${NODEPS}, $(findstring ${word}, ${MAKECMDGOALS}))))
-include $(patsubst %.f, %.dep, $(patsubst %.f90, %.dep, $(SRC_DEP)))
endif

# Build binaries for all physics models
allmodels:
	$(MAKE) MODEL=model199 clean
	$(MAKE) MODEL=model199 $(filter-out allmodels, ${MAKECMDGOALS})
#	$(MAKE) MODEL=model300 clean
#	$(MAKE) MODEL=model300 $(filter-out allmodels, ${MAKECMDGOALS})
#	$(MAKE) MODEL=model301 clean
#	$(MAKE) MODEL=model301 $(filter-out allmodels, ${MAKECMDGOALS})
	$(MAKE) MODEL=model302 clean
	$(MAKE) MODEL=model302 $(filter-out allmodels, ${MAKECMDGOALS})
#	$(MAKE) MODEL=model400 clean
#	$(MAKE) MODEL=model400 $(filter-out allmodels, ${MAKECMDGOALS})
#	$(MAKE) MODEL=model701 clean
#	$(MAKE) MODEL=model701 $(filter-out allmodels, ${MAKECMDGOALS})

#######################################################################
# Forcheck
#######################################################################
COMMA =,
SPACE :=
SPACE +=
OMPI_SRC = /afs/@cell/common/soft/intel/Compiler/11.1/c/include/omp_lib.f
OMPI_INC = $(dir $(OMPI_SRC))

# to remove the definitions 
# and evaluate the the shell scripts that return the parameters for pastix
INC_TMP = $(shell echo $(patsubst -D%,,$(INCLUDES))) 

DEFINES_FCK = $(DEFINES) -DFORCHECK

FCK_INCS = $(strip $(subst -I,,$(subst $(SPACE)-I,$(COMMA),$(INC_TMP))))
FCK_DEFS = $(strip $(subst -D,,$(subst $(SPACE)-D,$(COMMA),$(strip $(DEFINES_FCK)))))
FCK_SRC  = $(patsubst %.c,,$(JOREK2_MAIN_SRC))
FCK_POINCARE_SRC  = $(patsubst %.c,,$(JOREK2_POINCARE_SRC))
FCK_RST_BIN2HDF5_SRC = $(patsubst %.c,,$(RST_BIN2HDF5_SRC))
FCK_RST_HDF52BIN_SRC = $(patsubst %.c,,$(RST_HDF52BIN_SRC))
FCK_FOUR_SRC  = $(patsubst %.c,,$(JOREK2_FOUR_SRC))
FCK_POSTPROC_SRC  = $(patsubst %.c,,$(JOREK2_POSTPROC_SRC))
FCK_CONNECTION2_SRC  = $(patsubst %.c,,$(JOREK2_CONNECTION2_SRC))
FCK_STRIKES_SRC  = $(patsubst %.c,,$(JOREK2_STRIKES_SRC))
FCK_JOREK2VTK_SRC  = $(patsubst %.c,,$(JOREK2VTK_SRC))
FCK_JOREK2FLVTK_SRC  = $(patsubst %.c,,$(JOREK2FLVTK3D_SRC))
FCK_JOREK2VTK3D_SRC  = $(patsubst %.c,,$(JOREK2VTK3D_SRC))
FCK_DIAGNO_SRC  = $(patsubst %.c,,$(JOREK2_DIAGNO_SRC))
FCK_JOREK2_IMPORT_PERTURBATION_SRC  = $(patsubst %.c,,$(JOREK2_IMPORT_PERTURBATION_SRC))
FCK_CALL = forchk -allc -ancmpl -anref -declare -dp -l jorek.lst -define $(FCK_DEFS) \
	-I $(FCK_INCS)

forcheck : version 
	$(FCK_CALL) $(FCK_SRC) $(FCKDIR)/share/forcheck/MPI.flb

forcheck_poincare : version 
	$(FCK_CALL) diagnostics/jorek2_poincare.f90 $(FCK_POINCARE_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_rst_bin2hdf5 : version
	$(FCK_CALL) diagnostics/rst_bin2hdf5.f90 $(FCK_RST_BIN2HDF5_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_rst_hdf52bin : version
	$(FCK_CALL) diagnostics/rst_hdf52bin.f90 $(FCK_RST_HDF52BIN_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_four : version
	$(FCK_CALL) diagnostics/jorek2_four.f90 $(FCK_FOUR_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_postproc : version
	$(FCK_CALL),$(OMPI_INC) \
	postproc/jorek2_postproc.f90 $(FCK_POSTPROC_SRC) $(OMPI_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_connection2 : version
	$(FCK_CALL) diagnostics/jorek2_connection2.f90 $(FCK_CONNECTION2_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_strikes : version
	$(FCK_CALL) diagnostics/jorek2_strikes_ordered.f90 $(FCK_STRIKES_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jordel : version
	$(FCK_CALL) diagnostics/jordel.f90  $(FCK_JORDEL_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorpol : version
	$(FCK_CALL) diagnostics/jorpol.f90  $(FCK_JORPOL_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2vtk : version
	$(FCK_CALL) diagnostics/jorek2vtk.f90  $(FCK_JOREK2VTK_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_fieldlines_vtk : version
	$(FCK_CALL) diagnostics/jorek2_fieldlines_vtk.f90  $(FCK_JOREK2FLVTK_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2vtk_3d : version
	$(FCK_CALL) diagnostics/jorek2vtk_3d.f90 $(FCK_JOREK2VTK3D_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_diagno : version
	$(FCK_CALL) diagnostics/jorek2_diagno.f90   $(FCK_DIAGNO_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek_to_helena : version 
	$(FCK_CALL) diagnostics/jorek_to_helena.f90 models/mod_constants.f90 \
	timing/trace.f90 communication/mpi_mod.f90 $(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2_powers : version
	$(FCK_CALL) diagnostics/jorek2_powers.f90  $(FCK_JOREK2_POWERS_SRC) \
        $(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2_import_perturbation : version
	$(FCK_CALL) diagnostics/jorek2_import_perturbation.f90  $(FCK_JOREK2_IMPORT_PERTURBATION_SRC) \
        $(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2_target2vtk : version
	$(FCK_CALL) diagnostics/jorek2_target2vtk.f90  $(FCK_JOREK2_TARGET2VTK_SRC) \
        $(FCKDIR)/share/forcheck/MPI.flb

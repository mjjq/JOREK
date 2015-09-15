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
            NEW_DIAG_DEMO PENNING_TEST JOREK2_PARTICLES

# Add the common sources to all these programs
$(foreach prog,$(PROGRAMS),$(eval $(prog)_SRC += $(ALL_BINARIES_SRC) $(PPPSRC)))

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
		| $(AWK) -v file="$(patsubst %.f90, %.o, $<)" '{print tolower($$2)".mod : "file}' >> $@.tmp || touch $@.tmp;	\
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
		| $(AWK) -v file="$(patsubst %.f, %.o, $<)" '{print tolower($$2)".mod : "file}' >> $@.tmp || touch $@.tmp;	\
	fi;
	-@$(SED) -e "s/murge.inc//g" -e "s/dmumps_struc.h//g" < $@.tmp > $@ || touch $@
	-@rm -f $@.tmp

$(MAIN) : version $(JOREK2_MAIN_OBJ)
	$(FC) $(FFLAGS_OMP) \
	$(JOREK2_MAIN_OBJ) \
	 -o $(MAIN) $(INCLUDES) $(LIBS)

jorek2_poincare : version diagnostics/jorek2_poincare.f90 $(JOREK2_POINCARE_OBJ) 
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_poincare.f90 -o diagnostics/jorek2_poincare.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_poincare.o $(JOREK2_POINCARE_OBJ) \
	 -o $(JOREK_DIR)/jorek2_poincare  $(LIBS)

rst_bin2hdf5 : version diagnostics/rst_bin2hdf5.f90 $(RST_BIN2HDF5_OBJ) 
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/rst_bin2hdf5.f90 -o diagnostics/rst_bin2hdf5.o
	$(FC) $(FFLAGS_OMP) diagnostics/rst_bin2hdf5.o $(RST_BIN2HDF5_OBJ) \
	-o $(JOREK_DIR)/rst_bin2hdf5 $(LIBS)

rst_hdf52bin : version diagnostics/rst_hdf52bin.f90 $(RST_HDF52BIN_OBJ) 
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/rst_hdf52bin.f90 -o diagnostics/rst_hdf52bin.o
	$(FC) $(FFLAGS_OMP) diagnostics/rst_hdf52bin.o $(RST_HDF52BIN_OBJ) \
	-o $(JOREK_DIR)/rst_hdf52bin $(LIBS)

jorek2_four : version diagnostics/jorek2_four.f90 $(JOREK2_FOUR_OBJ)
	$(FC) $(FFLAGS) $(INCLUDES) -c diagnostics/jorek2_four.f90 -o diagnostics/jorek2_four.o
	$(FC) $(FFLAGS_OMP) diagnostics/jorek2_four.o $(JOREK2_FOUR_OBJ) \
	 -o $(JOREK_DIR)/jorek2_four $(LIBS) $(LIBFFTW)

jorek_extract_data : version $(JOREK_EXTRACT_DATA_OBJ) diagnostics/jorek_extract_data.o
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2_stan.f90	\
	$(JOREK2VTK_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_stan $(INCLUDES) $(LIBS)

jorek2_fieldlines_vtk : diagnostics/jorek2_fieldlines_vtk.f90 $(JOREK2FLVTK_OBJ)
	$(FC) $(FFLAGS)                   		\
	diagnostics/jorek2_fieldlines_vtk.f90         	\
	$(JOREK2FLVTK_OBJ)				\
	 -o $(JOREK_DIR)/jorek2_fieldlines_vtk $(INCLUDES) $(LIBS)

jorek2vtk_3d : diagnostics/jorek2vtk_3d.f90 $(JOREK2VTK3D_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2vtk_3d.f90    \
	$(JOREK2VTK3D_OBJ)		\
	 -o $(JOREK_DIR)/jorek2vtk_3d $(INCLUDES) $(LIBS)

jorek2_diagno : diagnostics/jorek2_diagno.f90 $(JOREK2_DIAGNO_OBJ)
	$(FC) $(FFLAGS_NO_OMP)          \
	diagnostics/jorek2_diagno.f90   \
	$(JOREK2_DIAGNO_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_diagno $(INCLUDES) $(LIBS)

jorek_to_helena : diagnostics/jorek_to_helena.f90
	$(FC) diagnostics/jorek_to_helena.f90 -o jorek_to_helena 

import_eqdsk : util/import_eqdsk.f90
	$(FC) util/import_eqdsk.f90 -o import_eqdsk $(LIBS)

jorek2_target2vtk : diagnostics/jorek2_target2vtk.f90 $(JOREK2_TARGET2VTK_OBJ)
	$(FC) $(FFLAGS)                         \
        diagnostics/jorek2_target2vtk.f90       \
        $(JOREK2_TARGET2VTK_OBJ)                \
         -o $(JOREK_DIR)/jorek2_target2vtk $(INCLUDES) $(LIBS)

jorek2_powers : diagnostics/jorek2_powers.f90 $(JOREK2_POWERS_OBJ)
	$(FC) $(FFLAGS_NO_OMP)                  \
        diagnostics/jorek2_powers.f90           \
        $(JOREK2_POWERS_OBJ)      	        \
         -o $(JOREK_DIR)/jorek2_powers $(INCLUDES) $(LIBS)

jorek2_import_perturbation : jorek2_import_perturbation_new_flags
jorek2_import_perturbation_new_flags: FFLAGS += -DIMPORT_PERTURBATIONS
jorek2_import_perturbation_new_flags: jorek2_import_perturbation_tmp

jorek2_import_perturbation_tmp : diagnostics/jorek2_import_perturbation.f90 $(JOREK2_IMPORT_PERTURBATION_OBJ)
	$(FC) $(FFLAGS) 				\
        diagnostics/jorek2_import_perturbation.f90	\
        $(JOREK2_IMPORT_PERTURBATION_OBJ)		\
         -o $(JOREK_DIR)/jorek2_import_perturbation $(INCLUDES) $(LIBS)
         
jorek2_particles : particles/jorek2_particles.f90 $(JOREK2_PARTICLES_OBJ)
	$(FC) $(FFLAGS)                 \
	particles/jorek2_particles.f90 	\
	$(JOREK2_PARTICLES_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_particles $(INCLUDES) $(LIBS) $(LIBFFTW)

penning_test : non_regression_tests/standalone/penning_test/penning.f90 $(PENNING_TEST_OBJ)
	$(FC) $(FFLAGS) non_regression_tests/standalone/penning_test/penning.f90 \
	$(PENNING_TEST_OBJ) -o $(JOREK_DIR)/penning_test \
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

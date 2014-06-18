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
	vacuum				\
	refinement			\
	postproc			\
	tools

LIBS = $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)

NODEPS = clean cleanall cleandep forcheck forcheck_poincare			\
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

JOREK2_MAIN_SRC         	+= $(ALL_BINARIES_SRC) $(PPPSRC) jorek2_main.f90
JOREK2_POINCARE_SRC     	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_POSTPROC_SRC     	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_POVRAY_SRC     		+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_CONNECTION2_SRC  	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_STRIKES_SRC		+= $(ALL_BINARIES_SRC) $(PPPSRC)
ENBIGGEN_SRC            	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JORDEL_SRC              	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JORPOL_SRC              	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2VTK_SRC           	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2FLVTK_SRC	        	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2VTK3D_SRC         	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_FOUR_SRC         	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_DIAGNO_SRC       	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_TARGET2VTK_SRC   	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_POWERS_SRC           	+= $(ALL_BINARIES_SRC) $(PPPSRC)
JOREK2_IMPORT_PERTURBATION_SRC	+= $(ALL_BINARIES_SRC) $(PPPSRC)

SRC_DEP = $(sort $(JOREK2_MAIN_SRC) $(JOREK2_POINCARE_SRC) $(JOREK2_CONNECTION2_SRC) $(JOREK2_STRIKES_SRC) $(JOREK2VTK_SRC) $(JOREK2FLVTK_SRC) \
          $(JOREK2VTK3D_SRC) $(JOREK2_FOUR_SRC) $(JOREK2_POSTPROC_SRC) $(JORDEL_SRC) $(JORPOL_SRC) $(ENBIGGEN_SRC)       \
	  $(JOREK2_TARGET2VTK_SRC) $(JOREK2_POWERS_SRC) $(JOREK2_IMPORT_PERTURBATION_SRC))

SRC_DEP := $(filter %.f90, $(SRC_DEP)) $(filter %.f, $(SRC_DEP)) diagnostics/hdf5_library.important

JOREK2_MAIN_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_MAIN_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_MAIN_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_MAIN_SRC)))

JOREK2_POINCARE_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_POINCARE_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_POINCARE_SRC)))	\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_POINCARE_SRC)))

JOREK2_CONNECTION2_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_CONNECTION2_SRC))) 	\
				$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_CONNECTION2_SRC)))		\
				$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_CONNECTION2_SRC)))

JOREK2_STRIKES_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_STRIKES_SRC))) 	\
				$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_STRIKES_SRC)))		\
				$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_STRIKES_SRC)))

ENBIGGEN_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(ENBIGGEN_SRC))) 	\
		$(patsubst %.f,%.o,$(filter %.f, $(ENBIGGEN_SRC)))      \
		$(patsubst %.c,%.o,$(filter %.c, $(ENBIGGEN_SRC)))

JORDEL_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JORDEL_SRC))) 	\
		$(patsubst %.f,%.o,$(filter %.f, $(JORDEL_SRC)))     	\
		$(patsubst %.c,%.o,$(filter %.c, $(JORDEL_SRC)))

JORPOL_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JORPOL_SRC))) 	\
		$(patsubst %.f,%.o,$(filter %.f, $(JORPOL_SRC)))     	\
		$(patsubst %.c,%.o,$(filter %.c, $(JORPOL_SRC)))

JOREK2VTK_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2VTK_SRC))) \
		$(patsubst %.f,%.o,$(filter %.f, $(JOREK2VTK_SRC)))	\
		$(patsubst %.c,%.o,$(filter %.c, $(JOREK2VTK_SRC)))

JOREK2FLVTK_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2FLVTK_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2FLVTK_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2FLVTK_SRC)))

JOREK2VTK3D_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2VTK3D_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2VTK3D_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2VTK3D_SRC)))

JOREK2_DIAGNO_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_DIAGNO_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_DIAGNO_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_DIAGNO_SRC)))

JOREK2_FOUR_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_FOUR_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_FOUR_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_FOUR_SRC)))

JOREK2_POSTPROC_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_POSTPROC_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_POSTPROC_SRC)))	\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_POSTPROC_SRC)))

JOREK2_POVRAY_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_POVRAY_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_POVRAY_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_POVRAY_SRC)))

JOREK2_TARGET2VTK_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_TARGET2VTK_SRC)))         \
                        $(patsubst %.f,%.o,$(filter %.f, $(JOREK2_TARGET2VTK_SRC)))             \
                        $(patsubst %.c,%.o,$(filter %.c, $(JOREK2_TARGET2VTK_SRC)))

JOREK2_POWERS_OBJ =     $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_POWERS_SRC)))     \
                        $(patsubst %.f,%.o,$(filter %.f, $(JOREK2_POWERS_SRC)))         \
                        $(patsubst %.c,%.o,$(filter %.c, $(JOREK2_POWERS_SRC)))

JOREK2_IMPORT_PERTURBATION_OBJ =	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_IMPORT_PERTURBATION_SRC)))     \
                        		$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_IMPORT_PERTURBATION_SRC)))         \
                        		$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_IMPORT_PERTURBATION_SRC)))

MOD_FILES=`find . -name "*.mod"`
MAIN = jorek_$(MODEL)

all: version $(MAIN)

cleanall : clean cleandep cleangenmod

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

version:
	@echo "#define SVN_VERSION"                                               > version.h.tmp
	@svn info > /dev/null 2>&1;                                                                \
	if [ $$? -eq 0 ]; then                                                                     \
	  export LANG=C; svn info | grep "Revision:" | sed -e 's/^Revision: *//' >> version.h.tmp; \
	else                                                                                       \
	  echo '"UNKNOWN"' >> version.h.tmp;                                                       \
	fi
	@cat version.h.tmp | tr '\n' ' ' > version.h
	@rm version.h.tmp
	@echo "" >> version.h
	@echo "#define compile_command '$(FC)'" >> version.h
	@echo "#define compile_flags '$(FFLAGS)'" >> version.h
	@echo "#define compile_includes '$(INCLUDES)'" >> version.h
	@echo "#define compile_defines '$(DEFINES)'" >> version.h
	@echo "#define compile_libs '$(LIBS)'" >> version.h
	-@echo "#define compile_dir '`pwd`'" >> version.h
	-@echo "#define compile_time '`date \"+%F %T\"`'" >> version.h
	-@echo "#define compile_user '`whoami`'" >> version.h
	-@echo "#define compile_machine '`hostname`'" >> version.h

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

$(MAIN) : $(JOREK2_MAIN_OBJ)
	$(FC) $(FFLAGS_OMP)	\
	$(JOREK2_MAIN_OBJ) 	\
	 -o $(MAIN) $(INCLUDES) $(LIBS)

jorek2_poincare : diagnostics/jorek2_poincare.f90 $(JOREK2_POINCARE_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2_poincare.f90 \
	$(JOREK2_POINCARE_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_poincare $(INCLUDES) $(LIBS)

jorek2_four : diagnostics/jorek2_four.f90 $(JOREK2_FOUR_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2_four.f90 	\
	$(JOREK2_FOUR_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_four $(INCLUDES) $(LIBS) $(LIBFFTW)

jorek2_postproc : postproc/jorek2_postproc.f90 $(JOREK2_POSTPROC_OBJ)
	$(FC) $(FFLAGS)                 \
	postproc/jorek2_postproc.f90 	\
	$(JOREK2_POSTPROC_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_postproc $(INCLUDES) $(LIBS) $(LIBFFTW)

jorek2_povray : diagnostics/jorek2_povray.f90 $(JOREK2_POVRAY_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2_povray.f90 	\
	$(JOREK2_POVRAY_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_povray $(INCLUDES) $(LIBS) $(LIBFFTW)

jorek2_connection2 : 	diagnostics/jorek2_connection2.f90 $(JOREK2_CONNECTION2_OBJ)
	$(FC) $(FFLAGS_OMP)                  	\
	diagnostics/jorek2_connection2.f90   	\
	$(JOREK2_CONNECTION2_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_connection2 $(INCLUDES) $(LIBS)

jorek2_connection_stan : 	diagnostics/jorek2_connection_stan.f90 $(JOREK2_CONNECTION2_OBJ)
	$(FC) $(FFLAGS_OMP)                  	\
	diagnostics/jorek2_connection_stan.f90   	\
	$(JOREK2_CONNECTION2_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_connection_stan $(INCLUDES) $(LIBS)

jorek2_strikes : 	diagnostics/jorek2_strikes_ordered.f90 $(JOREK2_STRIKES_OBJ)
	$(FC) $(FFLAGS_OMP)                  \
	diagnostics/jorek2_strikes_ordered.f90   	\
	$(JOREK2_STRIKES_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_strikes $(INCLUDES) $(LIBS)

enbiggen : diagnostics/enbiggen.f90 $(ENBIGGEN_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/enbiggen.f90        \
	$(ENBIGGEN_OBJ)              	\
	-o $(JOREK_DIR)/enbiggen $(INCLUDES) $(LIBS)

jordel : diagnostics/jordel.f90 $(JORDEL_OBJ)
	$(FC) $(FFLAGS)  		\
	diagnostics/jordel.f90          \
	$(JORDEL_OBJ)    		\
	-o $(JOREK_DIR)/jordel $(INCLUDES) $(LIBS)

jorpol : diagnostics/jorpol.f90 $(JORPOL_OBJ)
	$(FC) $(FFLAGS)  		\
	diagnostics/jorpol.f90          \
	$(JORPOL_OBJ)    		\
	-o $(JOREK_DIR)/jorpol $(INCLUDES) $(LIBS)

jorek2vtk : diagnostics/jorek2vtk.f90 $(JOREK2VTK_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2vtk.f90       \
	$(JOREK2VTK_OBJ)		\
	 -o $(JOREK_DIR)/jorek2vtk $(INCLUDES) $(LIBS)

jorek2_stan : diagnostics/jorek2_stan.f90 $(JOREK2VTK_OBJ)
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

forcheck :
	$(FCK_CALL) $(FCK_SRC) $(FCKDIR)/share/forcheck/MPI.flb

forcheck_poincare : 
	$(FCK_CALL) diagnostics/jorek2_poincare.f90 $(FCK_POINCARE_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_four :
	$(FCK_CALL) diagnostics/jorek2_four.f90 $(FCK_FOUR_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_postproc :
	$(FCK_CALL),$(OMPI_INC) \
	postproc/jorek2_postproc.f90 $(FCK_POSTPROC_SRC) $(OMPI_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_connection2 :
	$(FCK_CALL) diagnostics/jorek2_connection2.f90 $(FCK_CONNECTION2_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_strikes :
	$(FCK_CALL) diagnostics/jorek2_strikes_ordered.f90 $(FCK_STRIKES_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jordel :
	$(FCK_CALL) diagnostics/jordel.f90  $(FCK_JORDEL_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorpol :
	$(FCK_CALL) diagnostics/jorpol.f90  $(FCK_JORPOL_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2vtk :
	$(FCK_CALL) diagnostics/jorek2vtk.f90  $(FCK_JOREK2VTK_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_fieldlines_vtk :
	$(FCK_CALL) diagnostics/jorek2_fieldlines_vtk.f90  $(FCK_JOREK2FLVTK_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2vtk_3d :
	$(FCK_CALL) diagnostics/jorek2vtk_3d.f90 $(FCK_JOREK2VTK3D_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_diagno :
	$(FCK_CALL) diagnostics/jorek2_diagno.f90   $(FCK_DIAGNO_SRC) \
	$(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek_to_helena : 
	$(FCK_CALL) diagnostics/jorek_to_helena.f90 models/mod_constants.f90 \
	timing/trace.f90 communication/mpi_mod.f90 $(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2_powers :
	$(FCK_CALL) diagnostics/jorek2_powers.f90  $(FCK_JOREK2_POWERS_SRC) \
        $(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2_import_perturbation :
	$(FCK_CALL) diagnostics/jorek2_import_perturbation.f90  $(FCK_JOREK2_IMPORT_PERTURBATION_SRC) \
        $(FCKDIR)/share/forcheck/MPI.flb

forcheck_jorek2_target2vtk :
	$(FCK_CALL) diagnostics/jorek2_target2vtk.f90  $(FCK_JOREK2_TARGET2VTK_SRC) \
        $(FCKDIR)/share/forcheck/MPI.flb

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
	matrix				\
	solvers 			\
	plots				\
	diagnostics			\
	vacuum				\
	refinement			\
	tools

LIBS = $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)

NODEPS = clean cleanall cleandep
# If we have neither HIPS or PASTIX_MURGE we need to
# Use fake murge.

INCLUDES  =   -I. $(patsubst %,-I%/,$(DIRS)) $(INCMURGE)

VPATH = $(MAIN_MODEL_DIR) $(MODEL_DIR) $(DATATYPES_DIR) $(SOLVERS_DIR)

.SUFFIXES: .o .f90 .f .c


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

# Correct preprocessor-defines for IBM XLF Compiler
IBM_DEFINES = `echo $(DEFINES) | sed -e 's/^/-WF,/' -e 's/  */,/g'`

INCLUDES2  := $(INCLUDES) $(DEFINES)
ifdef IBMFC
  DEFINES  := $(DEFINES) -DIBM_MACHINE
  INCLUDES := $(INCLUDES) $(IBM_DEFINES)
else
  INCLUDES := $(INCLUDES) $(DEFINES)
endif

JOREK2_MAIN_SRC        = jorek2_main.f90 $(PPPSRC)
JOREK2_POINCARE_SRC    = $(PPPSRC)
JOREK2_CONNECTION2_SRC = $(PPPSRC)
JOREK2VTK_SRC          = $(PPPSRC)
JOREK2FLVTK_SRC	       = $(PPPSRC)
JOREK2VTK3D_SRC        = $(PPPSRC)
JOREK2_FOUR_SRC        = $(PPPSRC)

# include the description for each module
include $(patsubst %,%/module.mk,$(DIRS))

SRC_DEP = $(sort $(JOREK2_MAIN_SRC) $(JOREK2_POINCARE_SRC) $(JOREK2_CONNECTION2_SRC) $(JOREK2VTK_SRC) $(JOREK2FLVTK_SRC) $(JOREK2VTK3D_SRC) $(JOREK2_FOUR_SRC))

SRC_DEP := $(filter %.f90, $(SRC_DEP)) $(filter %.f, $(SRC_DEP))


JOREK2_MAIN_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_MAIN_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_MAIN_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_MAIN_SRC)))


JOREK2_POINCARE_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_POINCARE_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_POINCARE_SRC)))	\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_POINCARE_SRC)))

JOREK2_CONNECTION2_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_CONNECTION2_SRC))) 	\
				$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_CONNECTION2_SRC)))		\
				$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_CONNECTION2_SRC)))

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

MOD_FILES=`find . -name "*.mod"`
MAIN = jorek_$(MODEL)

all: version $(MAIN)

cleanall : clean cleandep

clean :	
	@echo ">> Deleting Object Files <<"
	-@rm -f $(JOREK2_MAIN_OBJ) $(JOREK2_FOUR_OBJ)
	@echo ">> Deleting Module Files <<"
	-@rm -f $(MOD_FILES)

cleandep:
	@echo ">> Deleting Dependency Files <<"
	-@rm -f *.dep */*.dep */*/*.dep

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
	-@echo "#define compile_time '`date \"+%F %T\"`'" >> version.h
	-@echo "#define compile_user '`whoami`'" >> version.h
	-@echo "#define compile_machine '`hostname`'" >> version.h

%.dep:%.f90
	@echo "Generating Dependencies for$(patsubst %.f90, %.o, $<)"
	@cpp $(INCLUDES2) < $< 2>/dev/null| grep -i "^[[:space:]]*use " | $(SED) 's/\,.*//' | 					\
		$(AWK) -v file_o=" $(patsubst %.f90, %.o, $<)" '{print file_o" : "tolower($$2)".mod"}' >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES2) < $< 2>/dev/null| 											\
		$(SED) -n "s@[ ]+include[ ]+'\([^']*\)*.*@$(patsubst %.f90, %.o, $<) : \1@pi" >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES2) < $< 2>/dev/null| 											\
		$(SED) -n 's@[ ]+include[ ]+"\([^"]*\)*.*@$(patsubst %.f90, %.o, $<) : \1@pi' >> $@.tmp || touch $@.tmp;
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
		| $(AWK) -v file="$(patsubst %.f, %.o, $<)" '{print tolower($$2)".mod : "file}' >> $@.tmp || touch $@.tmp;		\
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
	diagnostics/jorek2_four.f90 \
	$(JOREK2_FOUR_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_four $(INCLUDES) $(LIBS) $(LIBFFTW)

jorek2_connection2 : 	diagnostics/jorek2_connection2.f90 $(JOREK2_CONNECTION2_OBJ)
	$(FC) $(FFLAGS_OMP)                  \
	diagnostics/jorek2_connection2.f90   	\
	$(JOREK2_CONNECTION2_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_connection $(INCLUDES) $(LIBS)

jorek2vtk : diagnostics/jorek2vtk.f90 $(JOREK2VTK_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2vtk.f90       \
	$(JOREK2VTK_OBJ)		\
	 -o $(JOREK_DIR)/jorek2vtk $(INCLUDES) $(LIBS)

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
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2_diagno.f90   \
	$(JOREK2_DIAGNO_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_diagno $(INCLUDES) $(LIBS)

jorek_to_helena : diagnostics/jorek_to_helena.f90
	$(FC) diagnostics/jorek_to_helena.f90 -o jorek_to_helena 

import_eqdsk : util/import_eqdsk.f90
	$(FC) util/import_eqdsk.f90 -o import_eqdsk $(LIBS)

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

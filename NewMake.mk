include Makefile.inc

SED      ?=sed# gnu-sed command
AWK			 ?=awk# gnu-awk command
IBMFC    ?=#    IBM compiler flag (for FORTRAN symbol defs)


JOREK_DIR = `pwd`

DIRS =  datatypes models/$(MODEL) 	\
	models 				\
	communication 			\
	elements 			\
	grids 				\
	matrix 				\
	solvers 			\
	plots 				\
	diagnostics 			\
	vacuum 				\
	refinement			\
	timing				\
	tools

LIBS = $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)

NODEPS = clean cleanall cleandep
# If we have neither HIPS or PASTIX_MURGE we need to
# Use fake murge.

INCLUDES  =  $(patsubst %,-I%/,$(DIRS)) -I. $(INCMURGE)

VPATH = $(MAIN_MODEL_DIR) $(MODEL_DIR) $(DATATYPES_DIR) $(SOLVERS_DIR)

.SUFFIXES: .o .f90 .f .c


ifeq (1, $(USE_PASTIX_MURGE))
LIBS     := $(LIBS) $(LIB_PASTIX_MURGE) $(LIB_PASTIX_BLAS)
ifdef IBMFC 
	FDEFINES := $(DEFINES) -WF,-DUSE_MURGE
endif
DEFINES  := $(DEFINES) -DUSE_MURGE
INCLUDES := $(INCLUDES) $(INC_PASTIX)
endif

ifeq (1, $(USE_PASTIX))
DEFINES  := $(DEFINES) -DUSE_PASTIX
ifdef FDEFINES
	FDEFINES := $(FDEFINES),-DUSE_PASTIX
else 
	FDEFINES := -WF,-DUSE_PASTIX
endif
ifeq (0, $(USE_PASTIX_MURGE))
LIBS     := $(LIBS) $(LIB_PASTIX) $(LIB_PASTIX_BLAS)
INCLUDES := $(INCLUDES) $(INC_PASTIX)
endif
endif

ifeq (1, $(USE_HIPS))
LIBS := $(LIBS) $(LIBHIPS)
INCLUDES := $(INCLUDES) $(INCHIPS)
DEFINES  := $(DEFINES) -DUSE_HIPS
ifdef FDEFINES
	FDEFINES := $(FDEFINES),-DUSE_HIPS
else 
	FDEFINES := -WF,-DUSE_HIPS
endif
endif
ifeq (1, $(USE_MUMPS))
LIBS := $(LIBS) $(LIB_MUMPS) $(ORDLIB) $(SCALAP) $(BLACS) $(LIBLAPACK) $(LIBBLAS) $(PPPLIB) $(OPENMP_LIB)
INCLUDES := $(INCLUDES) -I$(INC_MUMPS)
DEFINES := $(DEFINES) -DUSE_MUMPS
ifdef FDEFINES
	FDEFINES := $(FDEFINES),-DUSE_MUMPS
else 
	FDEFINES := -WF,-DUSE_MUMPS
endif
endif

CINCLUDES := $(INCLUDES) $(DEFINES)
ifdef IBMFC
	INCLUDES  := $(INCLUDES) $(FDEFINES)
else
	INCLUDES  := $(CINCLUDES)
endif

JOREK2_MAIN_SRC        = jorek2_main.f90 $(PPPSRC)
JOREK2_POINCARE_SRC    = $(PPPSRC)
JOREK2_CONNECTION2_SRC = $(PPPSRC)
JOREK2VTK_SRC          = $(PPPSRC)
JOREK2FLVTK_SRC	       = $(PPPSRC)
JOREK2VTK3D_SRC        = $(PPPSRC)

# include the description for
#   each module
include $(patsubst %,%/module.mk,$(DIRS))

SRC_DEP = $(sort $(JOREK2_MAIN_SRC) $(JOREK2_POINCARE_SRC) $(JOREK2_CONNECTION2_SRC) $(JOREK2VTK_SRC) $(JOREK2FLVTK_SRC) $(JOREK2VTK3D_SRC))
SRC_DEP := $(filter %.f90, $(SRC_DEP)) $(filter %.f, $(SRC_DEP))


JOREK2_MAIN_OBJ = 	$(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_MAIN_SRC))) 	\
			$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_MAIN_SRC)))		\
			$(patsubst %.c,%.o,$(filter %.c, $(JOREK2_MAIN_SRC)))
JOREK2_POINCARE_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_POINCARE_SRC))) \
	$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_POINCARE_SRC)))
JOREK2_CONNECTION2_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_CONNECTION2_SRC))) \
	$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_CONNECTION2_SRC)))
JOREK2VTK_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2VTK_SRC))) \
	$(patsubst %.f,%.o,$(filter %.f, $(JOREK2VTK_SRC)))
JOREK2FLVTK_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2FLVTK_SRC))) \
	$(patsubst %.f,%.o,$(filter %.f, $(JOREK2FLVTK_SRC)))
JOREK2VTK3D_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2VTK3D_SRC))) \
	$(patsubst %.f,%.o,$(filter %.f, $(JOREK2VTK3D_SRC)))


MOD_FILES=`find . -name "*.mod"`
MAIN = jorek_$(MODEL)

all: $(MAIN)

modules :
	for dir in $(DIRS); do     \
          ($(MAKE) -C $$dir modules) \
        done

sources :	
	for dir in $(DIRS); do \
          ($(MAKE) -C $$dir all) \
        done

cleanall : clean cleandep

clean :	
	@echo ">> Deleting Object Files <<"
	-@rm -f $(JOREK2_MAIN_OBJ);
	@echo ">> Deleting Module Files <<"
	-@rm -f $(MOD_FILES);

cleandep:
	@echo ">> Deleting Dependency Files <<"
	-@rm -f */*.dep */*/*.dep;


%.dep:%.f90
	@echo "Generating Dependencies for $(patsubst %.f90, %.o, $<)"
	@cpp $(INCLUDES) < $< 2>/dev/null| grep -i "^[[:space:]]*use " | $(SED) 's/\,.*//' | 					\
		$(AWK) -v file_o=" $(patsubst %.f90, %.o, $<)" '{print file_o" : "tolower($$2)".mod"}' >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES) < $< 2>/dev/null| 											\
		$(SED) -n "s@[ ]+include[ ]+'\([^']*\)*.*@$(patsubst %.f90, %.o, $<) : \1@pi" >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES) < $< 2>/dev/null| 											\
		$(SED) -n 's@[ ]+include[ ]+"\([^"]*\)*.*@$(patsubst %.f90, %.o, $<) : \1@pi' >> $@.tmp || touch $@.tmp;
	@cpp $(INCLUDES) < $< 2>/dev/null|											\
		$(SED) -n 's@^\# *[0-9][0-9]* *"\([^"]*\)".*@$(patsubst %.f90, %.o, $<): \1@p' | 					\
		sort | uniq | grep -v ": /" | grep -v ": <">> $@.tmp || touch $@.tmp;
	@grep -q -i "^[[:space:]]*module" $< ; 											\
	if [ $$? -eq 0 ]; then 													\
		grep -i "^[[:space:]]*module" $<										\
		| $(AWK) -v file="$(patsubst %.f90, %.o, $<)" '{print tolower($$2)".mod : "file}' >> $@.tmp || touch $@.tmp;	\
	fi;
	-@$(SED) -e "s/murge.inc//g" -e "s/dmumps_struc.h//g" < $@.tmp > $@ || touch $@
	-@rm -f $@.tmp

%.dep: %.f
	@echo "Generating Dependencies for $(patsubst %.f, %.o, $<)"
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
	$(JOREK2_MAIN_OBJ)	\
	 -o $(MAIN) $(INCLUDES) $(LIBS)

jorek2_poincare : diagnostics/jorek2_poincare.f90 $(JOREK2_POINCARE_OBJ)
	$(FC) $(FFLAGS)                 \
	diagnostics/jorek2_poincare.f90 \
	$(JOREK2_POINCARE_OBJ)		\
	 -o $(JOREK_DIR)/jorek2_poincare $(INCLUDES) $(LIBS)

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
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(patsubst %.f, %.dep, $(patsubst %.f90, %.dep, $(SRC_DEP)))
endif

# Build binaries for all physics models
allmodels:
	$(MAKE) -f NewMake.mk MODEL=model199 clean
	$(MAKE) -f NewMake.mk MODEL=model199 $(filter-out allmodels, ${MAKECMDGOALS})
	
#	$(MAKE) -f NewMake.mk MODEL=model300 clean
#	$(MAKE) -f NewMake.mk MODEL=model300 $(filter-out allmodels, ${MAKECMDGOALS})
	
#	$(MAKE) -f NewMake.mk MODEL=model301 clean
#	$(MAKE) -f NewMake.mk MODEL=model301 $(filter-out allmodels, ${MAKECMDGOALS})
	
	$(MAKE) -f NewMake.mk MODEL=model302 clean
	$(MAKE) -f NewMake.mk MODEL=model302 $(filter-out allmodels, ${MAKECMDGOALS})
	
#	$(MAKE) -f NewMake.mk MODEL=model400 clean
#	$(MAKE) -f NewMake.mk MODEL=model400 $(filter-out allmodels, ${MAKECMDGOALS})
	
#	$(MAKE) -f NewMake.mk MODEL=model701 clean
#	$(MAKE) -f NewMake.mk MODEL=model701 $(filter-out allmodels, ${MAKECMDGOALS})

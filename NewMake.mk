include Makefile.inc

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
	timing

LIBS = $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)

NODEPS = clean cleanall cleandep
# If we have neither HIPS or PASTIX_MURGE we need to
# Use fake murge.

INCLUDES  =  $(patsubst %,-I%/,$(DIRS)) -I. $(INCMURGE)

VPATH = $(MAIN_MODEL_DIR) $(MODEL_DIR) $(DATATYPES_DIR) $(SOLVERS_DIR)

.SUFFIXES: .o .f90 .f 


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
LIBS := $(LIBS) $(LIB_MUMPS) $(ORDLIB) $(SCALAP) $(BLACS) $(LIBLAPACK) $(LIBBLAS) $(PPPLIB) $(OPENMP_LIB)
INCLUDES := $(INCLUDES) -I$(INC_MUMPS)
DEFINES := $(DEFINES) -DUSE_MUMPS
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

SRC_DEP = $(JOREK2_MAIN_SRC) $(JOREK2_POINCARE_SRC) $(JOREK2_CONNECTION2_SRC) $(JOREK2VTK_SRC) $(JOREK2FLVTK_SRC) $(JOREK2VTK3D_SRC)
SRC_DEP := $(shell echo "$(SRC_DEP)" | sed -e 's@ @\n@g' | sort -u)
SRC_DEP := $(filter %.f90, $(SRC_DEP)) $(filter %.f, $(SRC_DEP))


JOREK2_MAIN_OBJ = $(patsubst %.f90,%.o,$(filter %.f90, $(JOREK2_MAIN_SRC))) \
	$(patsubst %.f,%.o,$(filter %.f, $(JOREK2_MAIN_SRC)))
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


MOD_FILES=$(shell grep -i "^[[:space:]]*module" $(JOREK2_MAIN_SRC) | awk '{print tolower($$2)".mod"}' | xargs echo)
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
	@echo ">> suppression des objets pour jorek2_main <<"
	-@rm $(JOREK2_MAIN_OBJ);
	@echo ">> suppression des modules pour jorek2_main <<"
	-@rm $(MOD_FILES);

cleandep:
	echo ">> suppression des dependances pour jorek2_main <<"
	-@rm $(JOREK2_MAIN_DEP);

dependencies.mk: $(SRC_DEP)
	for file in $(SRC_DEP); do 				\
		file_o=`echo $$file|sed -e "s@\.f90@\.o@" |sed -e "s@\.f@\.o@" `; 	\
		echo $$file_o;							\
		cpp $(INCLUDES) < $$file | 						\
		grep -i "^[[:space:]]*use " | 						\
		awk -v file_o=$$file_o '{print file_o" : "tolower($$2)".mod"}' >> dependencies.tmp;	\
		cpp $(INCLUDES) < $$file | 						\
		sed -n "s@[ ]+include[ ]+'\([^']*\)*.*@"$$file_o" : \1@pi" >> dependencies.tmp; \
		cpp $(INCLUDES) < $$file | 						\
		sed -n 's@[ ]+include[ ]+"\([^"]*\)*.*@'$$file_o' : \1@pi' >> dependencies.tmp; \
		cpp $(INCLUDES) < $$file |						\
		 sed -n 's@^\# *[0-9][0-9]* *"\([^"]*\)".*@'$$file_o': \1@p' | \
		sort | uniq | grep -v ": /" | grep -v ": <">> dependencies.tmp; \
		file_o=`echo $$file | sed -e 's@f90@o@g'`;	\
		grep -q -i "^[[:space:]]*module" $$file ; 	\
		if [ $$? -eq 0 ]; then 				\
			grep -i "^[[:space:]]*module" $$file	\
			| awk -v file=$$file_o '{print tolower($$2)".mod : "file}' >> dependencies.tmp; \
		fi;						\
	done
	sed -e "s/murge.inc//g" -e "s/dmumps_struc.h//g" < dependencies.tmp > dependencies.mk
	rm dependencies.tmp

$(MAIN) : $(JOREK2_MAIN_OBJ) dependencies.mk
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
include dependencies.mk	
endif

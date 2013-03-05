DIR = models/model701
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)		\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/element_matrix.f90 		\
	$(DIR)/boundary_matrix_open.f90  	\
	$(DIR)/boundary_conditions.f90  	\
	$(DIR)/initial_conditions.f90           \
	$(DIR)/init_live_data_model.f90         \
	$(DIR)/mod_parameters.f90		

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/initialise_parameters.f90 	\
	$(DIR)/mod_parameters.f90

JOREK2_TARGET2VTK_SRC := $(JOREK2_TARGET2VTK_SRC) \
        $(DIR)/initialise_parameters.f90        \
        $(DIR)/mod_parameters.f90

JOREK2_POWERS_SRC := $(JOREK2_POWERS_SRC) \
        $(DIR)/initialise_parameters.f90        \
        $(DIR)/mod_parameters.f90

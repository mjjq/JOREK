DIR = models/model301

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/initialise_parameters.f90		\
	$(DIR)/mod_parameters.f90			\
	$(DIR)/sources.f90

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)		\
	$(DIR)/element_matrix.f90 		\
	$(DIR)/element_matrix_fft.f90 		\
	$(DIR)/boundary_matrix_open.f90  	\
	$(DIR)/boundary_matrix.f90 		\
	$(DIR)/boundary_conditions.f90 		\
	$(DIR)/initial_conditions.f90		\
	$(DIR)/init_live_data_model.f90

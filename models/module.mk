DIR=models
JOREK2_MAIN_SRC  := $(JOREK2_MAIN_SRC)			\
	$(DIR)/initialise_and_broadcast_parameters.f90	\
	$(DIR)/current.f90 				\
	$(DIR)/density.f90 				\
	$(DIR)/element_matrix_GS.f90 			\
	$(DIR)/element_matrix_GS_inverse.f90 		\
	$(DIR)/element_matrix_GS_perturbation.f90 	\
	$(DIR)/element_matrix_poisson.f90 		\
	$(DIR)/element_matrix_poisson_inverse.f90 	\
	$(DIR)/equilibrium.f90 				\
	$(DIR)/ffprime.f90 				\
	$(DIR)/poisson.f90 				\
	$(DIR)/temperature.f90 				\
	$(DIR)/temperature_i.f90 			\
	$(DIR)/temperature_e.f90 			\
	$(DIR)/mod_phys_module.f90 			\
	$(DIR)/mod_nodes_elements.f90			\
	$(DIR)/chgmt_node.f90				\
	$(DIR)/log_parameters.f90			\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/pellet_source.f90			\
	$(DIR)/derive_num_profiles.f90                  \
        $(DIR)/mod_pellet.f90


JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/derive_num_profiles.f90                  \
	$(DIR)/mod_phys_module.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/derive_num_profiles.f90                  \
	$(DIR)/mod_phys_module.f90


JOREK2VTK_SRC := $(JOREK2VTK_SRC)		\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/derive_num_profiles.f90                  \
	$(DIR)/mod_phys_module.f90

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)		\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/derive_num_profiles.f90                  \
	$(DIR)/mod_phys_module.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)		\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/derive_num_profiles.f90                  \
	$(DIR)/mod_phys_module.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/read_num_profiles.f90 			\
	$(DIR)/derive_num_profiles.f90                  \
	$(DIR)/mod_phys_module.f90                 \
	$(DIR)/density.f90 				\
	$(DIR)/temperature.f90 				\
        $(DIR)/mod_pellet.f90

DIR = models

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
        $(DIR)/current.f90 				\
	$(DIR)/density.f90 				\
	$(DIR)/derive_num_profiles.f90			\
	$(DIR)/det_modes.f90				\
	$(DIR)/ffprime.f90 				\
	$(DIR)/initialise_and_broadcast_parameters.f90	\
	$(DIR)/log_parameters.f90			\
	$(DIR)/mod_constants.f90			\
	$(DIR)/mod_domains.f90                          \
        $(DIR)/mod_equil_info.f90			\
	$(DIR)/mod_nodes_elements.f90			\
	$(DIR)/mod_phys_module.f90			\
	$(DIR)/preset_parameters.f90			\
	$(DIR)/read_num_profiles.f90			\
        $(DIR)/temperature.f90                          \
	$(DIR)/temperature_e.f90 			\
	$(DIR)/temperature_i.f90 			\
	$(DIR)/update_time_evol_params.f90

JOREK2_MAIN_SRC  := $(JOREK2_MAIN_SRC)			\
        $(DIR)/mod_bootstrap_functions.f90 			\
	$(DIR)/element_matrix_GS.f90 			\
	$(DIR)/element_matrix_GS_inverse.f90 		\
	$(DIR)/element_matrix_GS_perturbation.f90 	\
	$(DIR)/element_matrix_poisson.f90 		\
	$(DIR)/element_matrix_poisson_inverse.f90 	\
        $(DIR)/mod_corr_neg.f90                         \
	$(DIR)/equilibrium.f90 				\
	$(DIR)/F_profile.f90				\
	$(DIR)/poisson.f90 				\
	$(DIR)/chgmt_node.f90				\
	$(DIR)/pellet_source.f90			\
        $(DIR)/mod_diffusivities.f90                    \
        $(DIR)/mod_pellet.f90				\
        $(DIR)/read_RMP_profiles.f90                    \
        $(DIR)/neo_coef.f90				\
	$(DIR)/mod_assembly.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)			\
        $(DIR)/mod_pellet.f90				\

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
        $(DIR)/mod_pellet.f90				\

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)		\
        $(DIR)/mod_pellet.f90				\

JOREK2_POVRAY_SRC := $(JOREK2_POVRAY_SRC)		\
	$(DIR)/mod_pellet.f90 

RST_BIN2HDF5_SRC := $(RST_BIN2HDF5_SRC)                 \
	$(DIR)/mod_pellet.f90

RST_HDF52BIN_SRC := $(RST_HDF52BIN_SRC)                 \
	$(DIR)/mod_pellet.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
        $(DIR)/mod_pellet.f90             		\

JOREK2_STRIKES_SRC := $(JOREK2_STRIKES_SRC)      	\
	$(DIR)/mod_pellet.f90 

ENBIGGEN_SRC := $(JENBIGGEN_SRC)                        \
	$(DIR)/mod_pellet.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)			\
        $(DIR)/mod_corr_neg.f90                         \
        $(DIR)/mod_bootstrap_functions.f90 		\
	$(DIR)/mod_diffusivities.f90                    \
        $(DIR)/mod_pellet.f90                           \
        $(DIR)/neo_coef.f90                             

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)			\
        $(DIR)/mod_pellet.f90                           \

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)			\
        $(DIR)/mod_pellet.f90                           \

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)		\
	$(DIR)/mod_pellet.f90

JORDEL_SRC := $(JORDEL_SRC)                       	\
	$(DIR)/mod_diffusivities.f90                    \
        $(DIR)/mod_pellet.f90				\
	$(DIR)/neo_coef.f90

JORPOL_SRC := $(JORPOL_SRC)                       	\
	$(DIR)/mod_diffusivities.f90                    \
	$(DIR)/mod_pellet.f90                           \
	$(DIR)/neo_coef.f90

JOREK2_TARGET2VTK_SRC := $(JOREK2_TARGET2VTK_SRC)       \
        $(DIR)/mod_diffusivities.f90                    \
        $(DIR)/mod_pellet.f90 

JOREK2_POWERS_SRC := $(JOREK2_POWERS_SRC)               \
        $(DIR)/mod_diffusivities.f90                    \
        $(DIR)/mod_pellet.f90

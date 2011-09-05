DIR      = communication
JOREK2_MAIN_SRC  := $(JOREK2_MAIN_SRC)		\
	$(DIR)/broadcast_elements.f90 		\
	$(DIR)/broadcast_boundary.f90 		\
	$(DIR)/broadcast_nodes.f90 		\
	$(DIR)/broadcast_phys.f90 		\
	$(DIR)/broadcast_num_profiles.f90	\
	$(DIR)/distribute_harmonics.f90 	\
	$(DIR)/distribute_nodes_elements.f90 	\
	$(DIR)/distribute_vector.f90  		\
	$(DIR)/export_nemec.f90 		\
	$(DIR)/export_grid.f90 			\
	$(DIR)/export_pov.f90 			\
	$(DIR)/export_restart.f90 		\
	$(DIR)/import_restart.f90 		\
	$(DIR)/update_deltas.f90 		\
	$(DIR)/update_values.f90 		\
	$(DIR)/export_boundary.f90 		\
	$(DIR)/vertex_is_local.f90	 	\
	$(DIR)/export_helena.f90		\
	$(DIR)/split_broadcast.f90		\
        $(DIR)/mod_live_data.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)	\
	$(DIR)/broadcast_elements.f90 		\
	$(DIR)/broadcast_boundary.f90 		\
	$(DIR)/broadcast_nodes.f90 		\
	$(DIR)/broadcast_phys.f90 		\
	$(DIR)/export_boundary.f90 		\
	$(DIR)/import_restart.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/import_restart.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/import_restart.f90


JOREK2VTK_SRC := $(JOREK2VTK_SRC)		\
	$(DIR)/import_restart.f90

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)		\
	$(DIR)/import_restart.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)		\
	$(DIR)/import_restart.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/import_restart.f90		\
	$(DIR)/export_helena.f90

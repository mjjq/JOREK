DIR=vacuum
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 		\
	$(DIR)/ideal_wall.f90 			\
	$(DIR)/ideal_wall_starwall.f90		\
	$(DIR)/resistive_wall_starwall.f90	\
	$(DIR)/vacuum_old.f90			\
	$(DIR)/mod_vacuum.f90			\
	$(DIR)/mod_vacuum_response.f90		\
	$(DIR)/mod_vacuum_equilibrium.f90		

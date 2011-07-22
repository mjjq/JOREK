DIR=vacuum
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 		\
	$(DIR)/ideal_wall.f90 			\
	$(DIR)/ideal_wall_starwall.f90		\
	$(DIR)/resistive_wall_starwall.f90	\
	$(DIR)/vacuum_old.f90			\
	$(DIR)/mod_vacuum.f90			\
	$(DIR)/mod_vacuum_response.f90		\
	$(DIR)/mod_vacuum_equilibrium.f90		

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/mod_vacuum.f90			\

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/mod_vacuum.f90			\

JOREK2VTK_SRC := $(JOREK2VTK_SRC)		\
	$(DIR)/mod_vacuum.f90			\

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)		\
	$(DIR)/mod_vacuum.f90			\

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)		\
	$(DIR)/mod_vacuum.f90			\

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/mod_vacuum.f90			\

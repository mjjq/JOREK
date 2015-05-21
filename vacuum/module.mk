DIR = vacuum

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/mod_vacuum.f90

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 			\
	$(DIR)/mod_vacuum_response.f90			\
	$(DIR)/mod_vacuum_equilibrium.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)			\
	$(DIR)/mod_vacuum_response.f90			\
	$(DIR)/mod_vacuum_equilibrium.f90		

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
	$(DIR)/mod_vacuum_response.f90			\
	$(DIR)/mod_vacuum_equilibrium.f90

NEW_DIAG_DEMO_SRC := $(NEW_DIAG_DEMO_SRC)		\
	$(DIR)/mod_vacuum_response.f90			\
	$(DIR)/mod_vacuum_equilibrium.f90

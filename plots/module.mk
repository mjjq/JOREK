DIR = plots

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/ppplib.f

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)		\
	$(DIR)/plot_flux_surfaces.f90   	\
	$(DIR)/plot_grid.f90            	\
	$(DIR)/plot_profiles.f90        	\
	$(DIR)/plot_solution.f90        	\
	$(DIR)/plot_velocity_profile.f90 	\
	$(DIR)/plot_coils.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/plot_velocity_profile.f90	\
	$(DIR)/plot_profiles.f90

DIR=grids
JOREK2_MAIN_SRC:=$(JOREK2_MAIN_SRC)	\
	$(DIR)/define_boundary.f90 	\
	$(DIR)/fft.f90 			\
	$(DIR)/fgauss.f90 		\
	$(DIR)/find_crossing.f90	\
	$(DIR)/find_R_surface.f90 	\
	$(DIR)/find_RZ.f90 		\
	$(DIR)/find_theta_surface.f90 	\
	$(DIR)/find_Z_surface.f90 	\
	$(DIR)/grid_bezier_square.f90 	\
	$(DIR)/grid_flux_surface.f90 	\
	$(DIR)/grid_polar_bezier.f90 	\
	$(DIR)/grid_polar_bezier_square.f90 \
	$(DIR)/grid_xpoint.f90 		\
	$(DIR)/meshac.f90 		\
	$(DIR)/spline_spwert.f90	\
	$(DIR)/spline_tb15a.f90 	\
	$(DIR)/spline_tg02a.f90		\
	$(DIR)/boundary_from_grid.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/find_RZ.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/find_RZ.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/find_RZ.f90			\
	$(DIR)/meshac.f90			\
	$(DIR)/fgauss.f90
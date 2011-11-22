DIR=grids
JOREK2_MAIN_SRC:=$(JOREK2_MAIN_SRC)	    \
	$(DIR)/create_new_node.f90 	    \
	$(DIR)/create_x_node.f90 	    \
	$(DIR)/define_boundary.f90 	    \
	$(DIR)/define_flux_values.f90 	    \
	$(DIR)/define_new_grid_points.f90   \
	$(DIR)/define_final_grid.f90        \
	$(DIR)/fft.f90 			    \
	$(DIR)/fgauss.f90 		    \
	$(DIR)/find_crossing.f90	    \
	$(DIR)/find_R_surface.f90 	    \
	$(DIR)/find_RZ.f90 		    \
	$(DIR)/find_theta_surface.f90 	    \
	$(DIR)/find_Z_surface.f90 	    \
        $(DIR)/find_strategic_points.f90    \
	$(DIR)/grid_bezier_square.f90 	    \
	$(DIR)/grid_flux_surface.f90 	    \
	$(DIR)/grid_polar_bezier.f90 	    \
	$(DIR)/grid_polar_bezier_square.f90 \
	$(DIR)/grid_xpoint.f90 		    \
	$(DIR)/grid_double_xpoint.f90 	    \
	$(DIR)/meshac.f90 		    \
	$(DIR)/spline_spwert.f90	    \
	$(DIR)/spline_tb15a.f90 	    \
	$(DIR)/spline_tg02a.f90		    \
	$(DIR)/mod_grid_xpoint_data.f90	    \
	$(DIR)/mod_boundary.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)	\
	$(DIR)/find_RZ.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC) \
	$(DIR)/find_RZ.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/find_RZ.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/find_RZ.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/find_RZ.f90			\
	$(DIR)/meshac.f90			\
	$(DIR)/fgauss.f90

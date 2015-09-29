DIR = grids

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/mod_boundary.f90

JOREK2_MAIN_SRC:=$(JOREK2_MAIN_SRC)	    \
	$(DIR)/define_boundary.f90 	    \
	$(DIR)/grid_bezier_square.f90 	    \
	$(DIR)/grid_bezier_square_polar.f90 \
	$(DIR)/grid_double_xpoint.f90 	    \
	$(DIR)/grid_flux_surface.f90 	    \
	$(DIR)/grid_polar_bezier.f90 	    \
	$(DIR)/grid_xpoint.f90 		    \
        $(DIR)/grid_xpoint_wall.f90         \
	$(DIR)/mod_grid_xpoint_data.f90

PENNING_TEST_SRC += $(DIR)/define_boundary.f90   \
		    $(DIR)/grid_bezier_square.f90\
		    $(DIR)/grid_polar_bezier.f90 \
		    $(DIR)/grid_bezier_square_polar.f90

DIR = grids/grid_utils

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/find_RZ.f90 				\
	$(DIR)/find_RZ_nearby.f90

JOREK2_MAIN_SRC:=$(JOREK2_MAIN_SRC)	    			\
	$(DIR)/check_point_is_inside_wall.f90 	    		\
	$(DIR)/check_element_boundary.f90	\
	$(DIR)/update_neighbours.f90	\
	$(DIR)/create_new_node.f90 	    			\
	$(DIR)/create_x_node.f90 	    			\
	$(DIR)/define_flux_values.f90 	    			\
	$(DIR)/define_new_grid_points.f90   			\
	$(DIR)/define_final_grid.f90        			\
	$(DIR)/fft.f90 			    			\
	$(DIR)/fgauss.f90 		    			\
	$(DIR)/find_crossing.f90	    			\
	$(DIR)/find_extrapolation_points_central_part.f90	\
	$(DIR)/find_R_surface.f90 	    			\
        $(DIR)/find_strategic_points.f90    			\
	$(DIR)/find_theta_surface.f90 	    			\
	$(DIR)/find_wall_crossing.f90	    			\
	$(DIR)/find_wall_crossings_with_flux_surface.f90	\
	$(DIR)/find_Z_surface.f90 	    			\
	$(DIR)/meshac.f90 		    			\
	$(DIR)/py_plots_grids.f90 		    		\
	$(DIR)/mod_high_resolution_wall.f90 			\
	$(DIR)/reorder_flux_surfaces.f90    			\
	$(DIR)/spline_spwert.f90	    			\
	$(DIR)/spline_tb15a.f90 	    			\
	$(DIR)/spline_tg02a.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)				\
	$(DIR)/check_point_is_inside_wall.f90 	    		\
	$(DIR)/define_flux_values.f90 	    			\
	$(DIR)/fgauss.f90 		    			\
	$(DIR)/find_crossing.f90	    			\
	$(DIR)/find_extrapolation_points_central_part.f90	\
        $(DIR)/find_strategic_points.f90    			\
	$(DIR)/find_theta_surface.f90				\
	$(DIR)/find_wall_crossings_with_flux_surface.f90	\
	$(DIR)/find_Z_surface.f90 	    			\
	$(DIR)/meshac.f90 		    			\
	$(DIR)/mod_high_resolution_wall.f90 			\
	$(DIR)/py_plots_grids.f90 		    		\
	$(DIR)/reorder_flux_surfaces.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)			\
	$(DIR)/meshac.f90					\
	$(DIR)/fgauss.f90

JORDEL_SRC := $(JORDEL_SRC)               			\
	$(DIR)/find_theta_surface.f90

JORPOL_SRC := $(JORPOL_SRC)               			\
	$(DIR)/find_theta_surface.f90

ENBIGGEN_SRC := $(ENBIGGEN_SRC)               			\
	$(DIR)/find_theta_surface.f90

JOREK2_PARTICLES_SRC := $(JOREK2_PARTICLES_SRC)                 \
	$(DIR)/update_neighbours.f90	                        \
	$(DIR)/check_element_boundary.f90                       \
	$(DIR)/fft.f90 			    			\
	$(DIR)/fgauss.f90 		    			

PENNING_TEST_SRC += $(DIR)/update_neighbours.f90                \
		    $(DIR)/check_element_boundary.f90           \
		    $(DIR)/fft.f90                              \
		    $(DIR)/fgauss.f90                           \
	$(DIR)/meshac.f90 		    			\
	$(DIR)/spline_spwert.f90	    			\
	$(DIR)/spline_tb15a.f90 	    			\
	$(DIR)/spline_tg02a.f90 				

SIMON_PARTICLE_TEST_SRC += $(DIR)/update_neighbours.f90 	\
		    $(DIR)/check_element_boundary.f90           \
		    $(DIR)/fft.f90                              \
		    $(DIR)/fgauss.f90                           \
	$(DIR)/meshac.f90 		    			\
	$(DIR)/spline_spwert.f90	    			\
	$(DIR)/spline_tb15a.f90 	    			\
	$(DIR)/spline_tg02a.f90 				

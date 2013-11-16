include ../../config.in	

all : check_point_is_inside_wall.o		\
      create_new_node.o 			\
      create_x_node.o				\
      define_flux_values.o			\
      define_new_grid_points.o			\
      define_final_grid.o			\
      fft.o					\
      fgauss.o					\
      find_crossing.o				\
      find_R_surface.o				\
      find_RZ.o 				\
      find_strategic_points.o			\
      find_theta_surface.o			\
      find_wall_crossing.o			\
      find_wall_crossing_with_flux_surface.o	\
      find_Z_surface.o				\
      meshac.o					\
      reorder_flux_surfaces.o			\
      spline_spwert.o				\
      spline_tb15a.o				\
      spline_tg02a.o
      
modules : mod_high_resolution_wall.o
    
clean :
	rm -f *.o *mod
    

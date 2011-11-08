include ../config.in	

all : create_new_node.o \
      create_x_node.o \
      define_boundary.o \
      define_flux_values.o \
      define_new_grid_points.o \
      define_final_grid.o \
      fft.o \
      fgauss.o \
      find_crossing.o \
      find_R_surface.o \
      find_RZ.o \
      find_theta_surface.o \
      find_Z_surface.o \
      find_strategic_points.o \
      grid_bezier_square.o \
      grid_flux_surface.o \
      grid_polar_bezier.o \
      grid_polar_bezier_square.o \
      grid_xpoint.o \
      grid_double_xpoint.o \
      meshac.o \
      spline_spwert.o \
      spline_tb15a.o \
      spline_tg02a.o
      
modules :  mod_boundary.o  mod_grid_xpoint_data.o    
    
clean :
	rm -f *.o *mod
    

include ../config.in	

all : define_boundary.o 		\
      grid_bezier_square.o		\
      grid_bezier_polar_square.o	\
      grid_double_xpoint.o		\
      grid_flux_surface.o		\
      grid_polar_bezier.o		\
      grid_xpoint.o			\
      grid_xpoint_wall.o
      
modules :  mod_boundary.o  mod_grid_xpoint_data.o    
    
clean :
	rm -f *.o *mod
    

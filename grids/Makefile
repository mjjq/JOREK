include ../config.in	

all : define_boundary.o \
      fft.o \
      fgauss.o \
      find_crossing.o \
      find_R_surface.o \
      find_RZ.o \
      find_theta_surface.o \
      find_Z_surface.o \
      grid_bezier_square.o \
      grid_flux_surface.o \
      grid_polar_bezier.o \
      grid_polar_bezier_square.o \
      grid_xpoint.o \
      grid_xpoint_wall.o \
      meshac.o \
      spline_spwert.o \
      spline_tb15a.o \
      spline_tg02a.o
      
modules :  mod_boundary.o     
    
clean :
	rm -f *.o *mod
    

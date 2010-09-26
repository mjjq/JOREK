include ../config.in	

all : boundary_check.o         \
      dlength.o                \
      energy.o                 \
      find_axis.o              \
      find_flux_surfaces.o     \
      find_xpoint.o            \
      flux_surface_add_line.o  \
      flux_surface_add_point.o \
      print_grid.o             \
      psi_minmax.o             \
      q_profile.o              \
      RZ_minmax.o              \
      integrals.o              \
      integrals_3D.o          
      
modules : mod_fourier.o

clean :
	rm -f *.o *.mod

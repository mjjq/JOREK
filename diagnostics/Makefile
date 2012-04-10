include ../config.in	

all : boundary_check.o         \
      dlength.o                \
      energy.o                 \
      find_axis.o              \
      find_flux_surfaces.o     \
      find_xpoint.o            \
      find_limiter.o           \
      flux_surface_add_line.o  \
      flux_surface_add_point.o \
      print_grid.o             \
      psi_minmax.o             \
      q_profile.o              \
      RZ_minmax.o              \
      integrals.o              \
      integral_current.o       \
      integrals_3D.o           \
      output_saving.o          \
      HDF5_io.o

modules : mod_fourier.o \

clean :
	rm -f *.o *.mod

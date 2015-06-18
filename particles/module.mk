DIR = particles

JOREK2_MAIN_SRC:=$(JOREK2_MAIN_SRC)	\
  $(DIR)/mod_particles.f90           \
  $(DIR)/openadas.f90               \
	$(DIR)/update_particles.f90     	\
	$(DIR)/initialise_particles.f90 	\
	$(DIR)/redistribute_particles.f90 \
	$(DIR)/project_particles.f90  \
  $(DIR)/particles_vtk.f90 

JOREK2_PARTICLES_SRC :=$(JOREK2_PARTICLES_SRC) \
  $(DIR)/mod_particles.f90           \
  $(DIR)/openadas.f90                \
	$(DIR)/update_particles.f90	       \
	$(DIR)/initialise_particles.f90        \
	$(DIR)/redistribute_particles.f90      \
	$(DIR)/project_particles.f90  \
	$(DIR)/particles_vtk.f90


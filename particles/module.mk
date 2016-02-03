DIR = particles

JOREK2_PARTICLES_SRC :=$(JOREK2_PARTICLES_SRC) \
  $(DIR)/mod_particles.f90           \
  $(DIR)/openadas.f90                \
  $(DIR)/update_particles.f90	       \
  $(DIR)/calc_EB.f90 		\
  $(DIR)/initialise_particles.f90        \
  $(DIR)/redistribute_particles.f90      \
  $(DIR)/project_particles.f90  \
  $(DIR)/particles_vtk.f90

PENNING_TEST_SRC += $(DIR)/mod_particles.f90 		\
  $(DIR)/update_particles.f90 		\
  $(DIR)/calc_EB.f90 		\
  $(DIR)/mod_coordinate_transforms.f90

SIMON_PARTICLE_TEST_SRC += $(DIR)/mod_particles.f90 \
  $(DIR)/update_particles.f90 		\
  $(DIR)/calc_EB.f90 		\
  $(DIR)/guiding_center_position.f90 		\
  $(DIR)/mod_coordinate_transforms.f90

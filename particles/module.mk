DIR = particles

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC) \
  $(DIR)/mod_particles.f90           \
  $(DIR)/openadas.f90                \
  $(DIR)/update_particles.f90	       \
  $(DIR)/calc_EB.f90 		\
  $(DIR)/initialise_particles.f90        \
  $(DIR)/redistribute_particles.f90      \
  $(DIR)/particles_vtk.f90

JOREK2_PARTICLES_SRC :=$(JOREK2_PARTICLES_SRC) \
  $(DIR)/mod_particles.f90           \
  $(DIR)/mod_import_export_particles.f90 \
  $(DIR)/openadas.f90                \
  $(DIR)/update_particles.f90	       \
  $(DIR)/calc_EB.f90 		\
  $(DIR)/initialise_particles.f90        \
  $(DIR)/redistribute_particles.f90      \
  $(DIR)/mod_project_particles.f90  \
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

PROJECT_PARTICLES_VTK_SRC := $(PROJECT_PARTICLES_VTK_SRC) \
  $(DIR)/mod_particles.f90           \
  $(DIR)/mod_project_particles.f90 \
  $(DIR)/mod_import_export_particles.f90

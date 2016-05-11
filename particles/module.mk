DIR = particles

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC) \
  $(DIR)/mod_particles.f90           \
  $(DIR)/particles_vtk.f90

JOREK2_PARTICLES_SRC :=$(JOREK2_PARTICLES_SRC) \
  $(DIR)/openadas.f90                \
  $(DIR)/mod_coronal.f90                \
  $(DIR)/calc_EB.f90 		\
  $(DIR)/mod_initialise_particles.f90        \
  $(DIR)/redistribute_particles.f90      \
  $(DIR)/mod_import_export_particles.f90 \
  $(DIR)/mod_import_restart_linear.f90 \
  $(DIR)/mod_project_particles.f90  \
  $(DIR)/update_particles.f90	       \
  $(DIR)/mod_ionisation_recombination.f90 \
  $(DIR)/mod_redistribute_particles.f90

PENNING_TEST_SRC := $(PENNING_TEST_SRC) \
  $(DIR)/mod_coordinate_transforms.f90

SIMON_PARTICLE_TEST_SRC := $(SIMON_PARTICLE_TEST_SRC) \
  $(DIR)/guiding_center_position.f90 		\
  $(DIR)/mod_coordinate_transforms.f90

PROJECT_PARTICLES_VTK_SRC := $(PROJECT_PARTICLES_VTK_SRC) \
  $(DIR)/mod_project_particles.f90 \
  $(DIR)/mod_import_export_particles.f90

COUNT_PARTICLES_VTK_SRC := $(COUNT_PARTICLES_VTK_SRC) \
  $(DIR)/mod_import_export_particles.f90

DUMP_PARTICLES_VTK_SRC := $(DUMP_PARTICLES_VTK_SRC) \
  $(DIR)/mod_import_export_particles.f90

PARTICLE_FLUX_COORDINATES_SRC := $(PARTICLE_FLUX_COORDINATES_SRC) \
  $(DIR)/mod_import_export_particles.f90

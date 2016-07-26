DIR = particles

JOREK2_PARTICLES_SRC :=$(JOREK2_PARTICLES_SRC) \
  $(DIR)/mod_particles.f90 		      \
  $(DIR)/util/mod_openadas.f90                \
  $(DIR)/util/mod_coronal.f90                \
  $(DIR)/fields/mod_fields.f90 		     \
  $(DIR)/fields/jorek/mod_jorek_fields.f90   \
  $(DIR)/initialisation/mod_initialise_particles.f90        \
  $(DIR)/input_output/mod_particle_io.f90 \
  $(DIR)/integration/jorek/mod_import_restart_linear.f90 \
  $(DIR)/diagnostics/mod_project_particles.f90  \
  $(DIR)/diagnostics/particles_vtk.f90  \
  $(DIR)/pushers/mod_pusher.f90 \
  $(DIR)/pushers/boris/mod_boris.f90 \
  $(DIR)/ion_rec/mod_ionisation_recombination.f90 \
  $(DIR)/extra/mod_redistribute_particles.f90  \
  $(DIR)/diagnostics/mod_particle_diagnostics.f90

PENNING_TEST_SRC := $(PENNING_TEST_SRC) \
  $(DIR)/util/mod_coordinate_transforms.f90

SIMON_PARTICLE_TEST_SRC := $(SIMON_PARTICLE_TEST_SRC) \
  $(DIR)/util/guiding_center_position.f90 		\
  $(DIR)/util/mod_coordinate_transforms.f90

PROJECT_PARTICLES_VTK_SRC := $(PROJECT_PARTICLES_VTK_SRC) \
  $(DIR)/diagnostics/mod_project_particles.f90 \
  $(DIR)/input_output/mod_particle_io.f90

COUNT_PARTICLES_VTK_SRC := $(COUNT_PARTICLES_VTK_SRC) \
  $(DIR)/input_output/mod_particle_io.f90

DUMP_PARTICLES_VTK_SRC := $(DUMP_PARTICLES_VTK_SRC) \
  $(DIR)/diagnostics/particles_vtk.f90  \
  $(DIR)/input_output/mod_particle_io.f90

PARTICLE_FLUX_COORDINATES_SRC := $(PARTICLE_FLUX_COORDINATES_SRC) \
  $(DIR)/input_output/mod_particle_io.f90 \
  $(DIR)/diagnostics/mod_particle_diagnostics.f90

PARTICLE_FLUX_COORDINATE_DIFFUSION_SRC := $(PARTICLE_FLUX_COORDINATE_DIFFUSION_SRC) \
  $(DIR)/input_output/mod_particle_io.f90 \
  $(DIR)/diagnostics/mod_particle_diagnostics.f90

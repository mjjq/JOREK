DIR = tools

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 	\
	$(DIR)/fortran_pthread.c

JOREK2VTK_SRC := $(JOREK2VTK_SRC) \
	$(DIR)/mod_vtk.f90

RST_BIN2HDF5_SRC := $(RST_BIN2HDF5_SRC)  \
	$(DIR)/cla.f90

RST_HDF52BIN_SRC := $(RST_HDF52BIN_SRC)  \
	$(DIR)/cla.f90

PROJECT_PARTICLES_VTK_SRC := $(PROJECT_PARTICLES_VTK_SRC) \
	$(DIR)/mod_vtk.f90

COUNT_PARTICLES_VTK_SRC := $(COUNT_PARTICLES_VTK_SRC) \
	$(DIR)/mod_vtk.f90

DUMP_PARTICLES_VTK_SRC := $(DUMP_PARTICLES_VTK_SRC) \
	$(DIR)/mod_vtk.f90

JOREK2_PARTICLES_SRC := $(JOREK2_PARTICLES_SRC) \
	$(DIR)/mod_rng.f90 \
	$(DIR)/pcg_basic.c \
	$(DIR)/mod_pcg32.f90 \
	$(DIR)/mod_pcg32_rng.f90 \
	$(DIR)/mod_sobseq.f90 \
	$(DIR)/mod_sobseq_rng.f90 \
	$(DIR)/mod_random_seed.f90 \
	$(DIR)/fgetpid.c \
	$(DIR)/mod_sampling.f90

DIR = elements

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/basis_functions3.f90     \
	$(DIR)/initialise_basis.f90			\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90				\
	$(DIR)/interp_PRZ.f90       \
	$(DIR)/interp_PRZ_delta.f90       \
	$(DIR)/interp3_RZ.f90	\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/mod_gauss.f90

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)			\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)			\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK_EXTRACT_DATA_SRC := $(JOREK_EXTRACT_DATA_SRC)	\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

NEW_DIAG_DEMO_SRC := $(NEW_DIAG_DEMO_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_POVRAY_SRC := $(JOREK2_POVRAY_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)		\
	$(DIR)/hermite_1d.f90

RST_BIN2HDF5_SRC := $(RST_BIN2HDF5_SRC)                 \
	$(DIR)/hermite_1d.f90

RST_HDF52BIN_SRC := $(RST_HDF52BIN_SRC)                 \
	$(DIR)/hermite_1d.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/hermite_1d.f90

JOREK2_STRIKES_SRC := $(JOREK2_STRIKES_SRC)		\
	$(DIR)/hermite_1d.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)			\
	$(DIR)/hermite_1d.f90

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC) 			\
	$(DIR)/hermite_1d.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)			\
	$(DIR)/hermite_1d.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/hermite_elements.f90

JORDEL_SRC := $(JORDEL_SRC)                    		\
	$(DIR)/hermite_1d.f90

JORPOL_SRC := $(JORPOL_SRC)                    		\
	$(DIR)/hermite_1d.f90

ENBIGGEN_SRC := $(ENBIGGEN_SRC)               		\
	$(DIR)/hermite_1d.f90

JOREK2_TARGET2VTK_SRC := $(JOREK2_TARGET2VTK_SRC)	\
	$(DIR)/hermite_1d.f90

JOREK2_POWERS_SRC := $(JOREK2_POWERS_SRC)               \
        $(DIR)/hermite_1d.f90

JOREK2_IMPORT_PERTURBATION_SRC := $(JOREK2_IMPORT_PERTURBATION_SRC)	\
        $(DIR)/hermite_1d.f90

JOREK2_PARTICLES_SRC := $(JOREK2_PARTICLES_SRC)         \
       $(DIR)/hermite_1d.f90

PENNING_TEST_SRC := $(PENNING_TEST_SRC) \
       $(DIR)/hermite_1d.f90

SIMON_PARTICLE_TEST_SRC := $(SIMON_PARTICLE_TEST_SRC) \
       $(DIR)/hermite_1d.f90

PROJECT_PARTICLES_VTK_SRC := $(PROJECT_PARTICLES_VTK_SRC) \
       $(DIR)/hermite_1d.f90

COUNT_PARTICLES_VTK_SRC := $(COUNT_PARTICLES_VTK_SRC) \
       $(DIR)/hermite_1d.f90

DUMP_PARTICLES_VTK_SRC := $(DUMP_PARTICLES_VTK_SRC) \
       $(DIR)/hermite_1d.f90

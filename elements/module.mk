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
	$(DIR)/mod_gauss.f90 				\
       $(DIR)/hermite_1d.f90

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)			\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)			\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK_EXTRACT_DATA_SRC := $(JOREK_EXTRACT_DATA_SRC)	\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

NEW_DIAG_DEMO_SRC := $(NEW_DIAG_DEMO_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_POVRAY_SRC := $(JOREK2_POVRAY_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)		\
	$(DIR)/bezier_1d.f90				\
	$(DIR)/hermite_elements.f90

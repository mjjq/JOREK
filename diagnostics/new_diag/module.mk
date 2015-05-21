DIR = diagnostics/new_diag

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
        $(DIR)/mod_new_diag.f90                         \
	$(DIR)/mod_position.f90				\
	$(DIR)/mod_expression.f90			\
	$(DIR)/mod_straight_field_line.f90              \
	$(DIR)/mod_four_filter.f90                      \
        $(DIR)/mod_diag_output.f90

NEW_DIAG_DEMO_SRC := $(NEW_DIAG_DEMO_SRC)		\
        $(DIR)/mod_new_diag.f90                         \
	$(DIR)/mod_position.f90				\
	$(DIR)/mod_expression.f90			\
	$(DIR)/mod_straight_field_line.f90              \
	$(DIR)/mod_four_filter.f90                      \
        $(DIR)/mod_diag_output.f90

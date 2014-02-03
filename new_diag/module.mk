DIR = new_diag

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
        $(DIR)/mod_new_diag.f90                         \
	$(DIR)/mod_position.f90				\
	$(DIR)/mod_expression.f90			\
	$(DIR)/mod_straight_field_line.f90              \
	$(DIR)/mod_four_filter.f90                      \
        $(DIR)/mod_diag_output.f90

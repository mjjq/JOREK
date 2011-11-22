DIR = postproc
JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)	  \
	$(DIR)/mod_settings.f90                   \
	$(DIR)/mod_convert_character.f90          \
	$(DIR)/mod_parse_commands.f90             \
	$(DIR)/mod_postproc_help.f90              \
	$(DIR)/mod_exec_commands.f90

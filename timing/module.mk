DIR = timing

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/clock.f90				\
	$(DIR)/flush_it.f90				\
	$(DIR)/pastix_getmem.c				\
	$(DIR)/trace.f90

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)			\
	$(DIR)/r3_info.f90				\
	$(DIR)/r3_ctlk.c				\
	$(DIR)/flushc.c					\

DIR = timing
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)	\
	$(DIR)/r3_info.f90		\
	$(DIR)/r3_ctlk.c		\
	$(DIR)/flushc.c			\
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/clock.f90		\
	$(DIR)/pastix_getmem.c

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/clock.f90		\
	$(DIR)/pastix_getmem.c


JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/clock.f90		\
	$(DIR)/pastix_getmem.c


JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/clock.f90		\
	$(DIR)/pastix_getmem.c


JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/pastix_getmem.c


JOREK2VTK_SRC := $(JOREK2VTK_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/clock.f90		\
	$(DIR)/pastix_getmem.c


JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/pastix_getmem.c


JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/clock.f90		\
	$(DIR)/pastix_getmem.c


JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC) \
	$(DIR)/flush_it.f90		\
        $(DIR)/clock.f90                \
	$(DIR)/trace.f90		\
	$(DIR)/pastix_getmem.c


JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC) \
	$(DIR)/flush_it.f90		\
	$(DIR)/trace.f90		\
	$(DIR)/pastix_getmem.c

JORDEL_SRC := $(JORDEL_SRC) \
	$(DIR)/flush_it.f90             \
	$(DIR)/trace.f90                \
	$(DIR)/clock.f90                \
	$(DIR)/pastix_getmem.c

ENBIGGEN_SRC := $(ENBIGGEN_SRC) \
	$(DIR)/flush_it.f90             \
	$(DIR)/trace.f90                \
	$(DIR)/clock.f90                \
	$(DIR)/pastix_getmem.c

JOREK2_TARGET2VTK_SRC := $(JOREK2_TARGET2VTK_SRC) \
        $(DIR)/flush_it.f90             \
        $(DIR)/trace.f90                \
        $(DIR)/pastix_getmem.c

JOREK2_POWERS_SRC := $(JOREK2_POWERS_SRC) \
        $(DIR)/flush_it.f90             \
        $(DIR)/trace.f90                \
        $(DIR)/pastix_getmem.c

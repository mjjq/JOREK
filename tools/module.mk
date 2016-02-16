DIR = tools

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 	\
	$(DIR)/fortran_pthread.c

JOREK2VTK_SRC := $(JOREK2VTK_SRC) \
	$(DIR)/mod_vtk.f90

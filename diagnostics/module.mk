DIR = diagnostics

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC) \
	$(DIR)/find_axis.f90	        \
	$(DIR)/find_limiter.f90	        \
	$(DIR)/find_strike.f90	        \
	$(DIR)/find_xpoint.f90	        \
	$(DIR)/integrals_3D.f90	        \
	$(DIR)/mod_diagnostics.f90      \
	$(DIR)/psi_minmax.f90	        \
	$(DIR)/RZ_minmax.f90            \
	$(DIR)/hdf5_io.f90            \
	$(DIR)/jorek2help.f90            
  

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)	  \
	$(DIR)/boundary_check.f90         \
	$(DIR)/dlength.f90                \
	$(DIR)/energy.f90                 \
	$(DIR)/temp.f90                   \
	$(DIR)/find_flux_surfaces.f90     \
	$(DIR)/flux_surface_add_line.f90  \
	$(DIR)/flux_surface_add_point.f90 \
	$(DIR)/print_grid.f90             \
	$(DIR)/integral_current.f90       \
	$(DIR)/q_profile.f90              \
	$(DIR)/determine_q_profile.f90    \
	$(DIR)/determine_PhiN.f90         \
	$(DIR)/integrals.f90              \
	$(DIR)/output_saving.f90        

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)	  \
	$(DIR)/boundary_check.f90         \
	$(DIR)/mod_fourier.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)	  \
	$(DIR)/boundary_check.f90         \
	$(DIR)/flux_surface_add_line.f90  \
	$(DIR)/flux_surface_add_point.f90 \
	$(DIR)/determine_q_profile.f90    \
	$(DIR)/find_flux_surfaces.f90     \
	$(DIR)/integrals.f90

NEW_DIAG_DEMO_SRC := $(NEW_DIAG_DEMO_SRC)	  \
	$(DIR)/boundary_check.f90         \
	$(DIR)/flux_surface_add_line.f90  \
	$(DIR)/flux_surface_add_point.f90 \
	$(DIR)/determine_q_profile.f90    \
	$(DIR)/find_flux_surfaces.f90     \
	$(DIR)/integrals.f90

JOREK2_POVRAY_SRC := $(JOREK2_POVRAY_SRC)	  \
	$(DIR)/find_flux_surfaces.f90     \
	$(DIR)/flux_surface_add_line.f90  	\
	$(DIR)/flux_surface_add_point.f90

JOREK2_STRIKES_SRC := $(JOREK2_STRIKES_SRC) 	\
	$(DIR)/divertor_desc.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)		\
	$(DIR)/find_flux_surfaces.f90  		\
	$(DIR)/flux_surface_add_line.f90  	\
	$(DIR)/determine_q_profile.f90    	\
	$(DIR)/flux_surface_add_point.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)           \
	$(DIR)/find_flux_surfaces.f90  		\
	$(DIR)/flux_surface_add_line.f90  	\
	$(DIR)/flux_surface_add_point.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/find_flux_surfaces.f90  		\
	$(DIR)/flux_surface_add_line.f90  	\
	$(DIR)/flux_surface_add_point.f90  	\
	$(DIR)/integrals.f90

JORDEL_SRC := $(JORDEL_SRC)                     \
	$(DIR)/find_flux_surfaces.f90           \
	$(DIR)/flux_surface_add_line.f90        \
	$(DIR)/flux_surface_add_point.f90

JORPOL_SRC := $(JORPOL_SRC)                     \
	$(DIR)/find_flux_surfaces.f90           \
	$(DIR)/flux_surface_add_line.f90        \
	$(DIR)/flux_surface_add_point.f90

JOREK2_POWERS_SRC := $(JOREK2_POWERS_SRC)   \
	$(DIR)/q_profile.f90                \
	$(DIR)/determine_q_profile.f90

include config.in

DIRS = datatypes models/$(MODEL) models communication elements grids matrix solvers timing plots diagnostics vacuum refinement tools

MAIN = jorek_$(MODEL)

.PHONY : modules sources $(MAIN) jorek2_diagno jorek2vtk jorek2vtk_3d import_eqdsk

all :   modules sources $(MAIN)

modules :
	@for dir in $(DIRS); do     \
          ($(MAKE) -C $$dir modules) || exit 1; \
        done

sources :	
	@for dir in $(DIRS); do \
          ($(MAKE) -C $$dir all) || exit 1; \
        done

clean :	
	@rm -f $(MAIN)
	@for dir in $(DIRS); do   \
          ($(MAKE) -C $$dir clean) \
        done


$(MAIN) : jorek2_main.f90
	$(FC) $(FFLAGS)      \
	jorek2_main.f90      \
	datatypes/*.o        \
	models/$(MODEL)/*.o  \
	models/*.o           \
	communication/*.o    \
	elements/*.o         \
	grids/*.o            \
	matrix/*.o           \
	solvers/*.o          \
        timing/r3_info.o     \
        timing/flushc.o      \
	refinement/*.o       \
	plots/*.o            \
        vacuum/*.o           \
	diagnostics/*.o      \
	 -o $(MAIN) $(INCLUDES) $(LIBS)

jorek2_four : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2_four.f90   \
	datatypes/*.o        \
	models/$(MODEL)/*.o  \
	models/*.o           \
	communication/*.o    \
	elements/*.o         \
	grids/*.o            \
	matrix/*.o           \
	solvers/*.o          \
        timing/*.o           \
	refinement/*.o       \
	plots/*.o            \
        vacuum/*.o           \
	diagnostics/*.o      \
	 -o $(JOREK_DIR)/jorek2_four $(INCLUDES) $(LIBS)

jorek2_poincare : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2_poincare.f90   \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	elements/mod_basis_at_gaussian.o  \
	models/$(MODEL)/mod_phys_module.o \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
	elements/basis_functions2.o       \
        elements/hermite_1d.o             \
        elements/interp.o                 \
	elements/interp_RZ.o              \
	diagnostics/find_axis.o           \
	diagnostics/find_xpoint.o         \
	diagnostics/RZ_minmax.o           \
	grids/find_RZ.o                   \
	solvers/root.o                    \
	solvers/mnewtax.o                 \
        refinement/neighbours.o           \
        refinement/Ref_Add_Elements.o\
        refinement/Ref_Add_Node.o\
        refinement/Ref_boundary_node.o\
        refinement/Ref_Check_Neighb_Stat.o\
        refinement/Ref_Find_Constrained_Node.o\
        refinement/Refine_Element.o\
        refinement/Refine_Elem_List.o\
        refinement/Ref_Update_Index.o\
        refinement/Ref_Active_Node.o\
	 -o $(JOREK_DIR)/jorek2_poincare $(INCLUDES) $(LIBS)

jorek2_connection2 : modules sources
	$(FC) $(FFLAGS_NO_OMP)                   \
	diagnostics/jorek2_connection2.f90   \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	elements/mod_basis_at_gaussian.o  \
	models/$(MODEL)/mod_phys_module.o \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
	elements/basis_functions2.o       \
        elements/hermite_1d.o             \
        elements/interp.o                 \
	elements/interp_RZ.o              \
	diagnostics/find_axis.o           \
	diagnostics/find_xpoint.o         \
	diagnostics/RZ_minmax.o           \
	grids/find_RZ.o                   \
	solvers/root.o                    \
	solvers/mnewtax.o                 \
        refinement/neighbours.o           \
	 -o $(JOREK_DIR)/jorek2_connection $(INCLUDES) $(LIBS)

jorek2vtk : modules sources
	$(FC) $(FFLAGS_NOOMP)             \
	diagnostics/jorek2vtk.f90         \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	models/$(MODEL)/mod_phys_module.o \
	elements/mod_basis_at_gaussian.o  \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
        elements/basis_functions1.o       \
	elements/basis_functions2.o       \
        elements/interp.o                 \
	elements/interp_RZ.o              \
	diagnostics/find_axis.o           \
	diagnostics/find_xpoint.o         \
	solvers/mnewtax.o                 \
	 -o $(JOREK_DIR)/jorek2vtk $(INCLUDES) $(LIBS)

jorek2_fieldlines_vtk : modules sources
	$(FC) $(FFLAGS_NOOMP)             \
	diagnostics/jorek2_fieldlines_vtk.f90         \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	models/$(MODEL)/mod_phys_module.o \
	elements/mod_basis_at_gaussian.o  \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
	elements/basis_functions2.o       \
        elements/interp.o                 \
	elements/interp_RZ.o              \
	diagnostics/find_axis.o           \
	diagnostics/find_xpoint.o         \
	solvers/mnewtax.o                 \
        refinement/neighbours.o           \
	 -o $(JOREK_DIR)/jorek2_fieldlines_vtk $(INCLUDES) $(LIBS)

jorek2vtk_3d : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2vtk_3d.f90      \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	elements/mod_gauss.o              \
	models/$(MODEL)/mod_phys_module.o \
	elements/mod_basis_at_gaussian.o  \
	communication/import_restart.o    \
	elements/initialise_basis.o       \
	elements/basis_functions.o        \
	elements/basis_functions1.o       \
	elements/basis_functions2.o       \
	elements/interp.o                 \
	elements/interp_RZ.o              \
	 -o $(JOREK_DIR)/jorek2vtk_3d $(INCLUDES) $(LIBS)

jorek2_diagno : modules sources
	$(FC) $(FFLAGS)                   \
	diagnostics/jorek2_diagno.f90     \
	datatypes/mod_parameters.o        \
	datatypes/mod_data_structure.o    \
	communication/export_helena.o     \
	plots/plot_velocity_profile.o     \
	plots/plot_profiles.o             \
	diagnostics/find_flux_surfaces.o  \
	diagnostics/flux_surface_add_line.o  \
	diagnostics/flux_surface_add_point.o  \
	diagnostics/find_axis.o           \
	diagnostics/find_xpoint.o         \
	diagnostics/RZ_minmax.o           \
	diagnostics/psi_minmax.o          \
	diagnostics/integrals.o           \
	diagnostics/integrals_3D.o        \
	solvers/mnewtax.o                 \
	solvers/solvP3.o                  \
	solvers/solve_M2.o                \
	solvers/root.o                    \
	elements/*.o                      \
	models/$(MODEL)/mod_phys_module.o \
	communication/import_restart.o    \
	grids/find_RZ.o                   \
	grids/meshac.o                    \
	grids/fgauss.o                    \
	 -o $(JOREK_DIR)/jorek2_diagno $(INCLUDES) $(LIBS)

jorek_to_helena :
	$(FC) diagnostics/jorek_to_helena.f90 -o jorek_to_helena 

import_eqdsk :
	$(FC) util/import_eqdsk.f90 -o import_eqdsk $(LIBS)

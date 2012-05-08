DIR = solvers
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)	\
	$(DIR)/gmres_driver.f90		\
	$(DIR)/gmres_matrix_vector.f90 	\
	$(DIR)/gmres_precondition.f90 	\
	$(DIR)/dPackgmres.f  		\
	$(DIR)/initialise_mumps.f90 	\
	$(DIR)/initialise_pastix.f90 	\
	$(DIR)/mnewtax.f90 		\
	$(DIR)/root.f90 		\
	$(DIR)/solve_M2.f90 		\
	$(DIR)/solve_matrix_n.f90 	\
	$(DIR)/solve_mumps_all.f90 	\
	$(DIR)/solve_pastix_all.f90 	\
	$(DIR)/solvP3.f90 		\
	$(DIR)/update_rhs_n.f90 	\
	$(DIR)/solve_murge_all.f90	\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90		\
	$(DIR)/mod_wsmp.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)	\
	$(DIR)/root.f90				\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\
	$(DIR)/mnewtax.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
	$(DIR)/root.f90				\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\
	$(DIR)/solve_M2.f90 		\
	$(DIR)/solvP3.f90 		\
	$(DIR)/mnewtax.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/root.f90				\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\
	$(DIR)/mnewtax.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/root.f90					\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\
	$(DIR)/mnewtax.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)		\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\
	$(DIR)/mnewtax.f90              \
	$(DIR)/mod_wsmp.f90           


JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)		\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\
	$(DIR)/mnewtax.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)   \
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)	\
	$(DIR)/mnewtax.f90	                \
	$(DIR)/solvP3.f90       	        \
	$(DIR)/solve_M2.f90             	\
	$(DIR)/root.f90                    	\
	$(DIR)/mod_mumps.f90 		\
	$(DIR)/mod_pastix.f90		\
	$(DIR)/mod_murge.f90 		\

DIR = elements
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC)			\
	$(DIR)/basis_functions2.f90 			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions.f90 			\
	$(DIR)/bezier_1d.f90 				\
	$(DIR)/hermite_1d.f90 				\
	$(DIR)/hermite_elements.f90 			\
	$(DIR)/initialise_basis.f90 			\
	$(DIR)/interp.f90 				\
	$(DIR)/interp_RZ.f90				\
	$(DIR)/mod_gauss.f90            		\
	$(DIR)/mod_basis_at_gaussian.f90

JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC)			\
	$(DIR)/basis_functions2.f90 			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions.f90 			\
	$(DIR)/bezier_1d.f90 				\
	$(DIR)/hermite_1d.f90 				\
	$(DIR)/hermite_elements.f90 			\
	$(DIR)/initialise_basis.f90 			\
	$(DIR)/interp.f90 				\
	$(DIR)/interp_RZ.f90				\
	$(DIR)/mod_gauss.f90            		\
	$(DIR)/mod_basis_at_gaussian.f90

JOREK2_POSTPROC_SRC := $(JOREK2_POSTPROC_SRC)		\
	$(DIR)/basis_functions2.f90 			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions.f90 			\
	$(DIR)/bezier_1d.f90 				\
	$(DIR)/hermite_1d.f90 				\
	$(DIR)/hermite_elements.f90 			\
	$(DIR)/initialise_basis.f90 			\
	$(DIR)/interp.f90 				\
	$(DIR)/interp_RZ.f90				\
	$(DIR)/mod_gauss.f90            		\
	$(DIR)/mod_basis_at_gaussian.f90

JOREK2_POVRAY_SRC := $(JOREK2_POVRAY_SRC)		\
	$(DIR)/basis_functions2.f90 			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions.f90 			\
	$(DIR)/bezier_1d.f90 				\
	$(DIR)/hermite_1d.f90 				\
	$(DIR)/hermite_elements.f90 			\
	$(DIR)/initialise_basis.f90 			\
	$(DIR)/interp.f90 				\
	$(DIR)/interp_RZ.f90				\
	$(DIR)/mod_gauss.f90            		\
	$(DIR)/mod_basis_at_gaussian.f90

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)		\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/initialise_basis.f90			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90


JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/initialise_basis.f90			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90


JOREK2_STRIKES_SRC := $(JOREK2_STRIKES_SRC)		\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/initialise_basis.f90			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/hermite_1d.f90				\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90


JOREK2VTK_SRC := $(JOREK2VTK_SRC)			\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/hermite_1d.f90  	                	\
	$(DIR)/initialise_basis.f90			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)			\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/initialise_basis.f90			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)			\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/hermite_1d.f90 				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/initialise_basis.f90			\
	$(DIR)/basis_functions.f90			\
	$(DIR)/basis_functions1.f90			\
	$(DIR)/basis_functions2.f90			\
	$(DIR)/interp.f90				\
	$(DIR)/interp_RZ.f90

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)		\
	$(DIR)/mod_gauss.f90				\
	$(DIR)/mod_basis_at_gaussian.f90		\
	$(DIR)/basis_functions2.f90 			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions.f90 			\
	$(DIR)/bezier_1d.f90 				\
	$(DIR)/hermite_1d.f90 				\
	$(DIR)/hermite_elements.f90 			\
	$(DIR)/initialise_basis.f90 			\
	$(DIR)/interp.f90 				\
	$(DIR)/interp_RZ.f90

JORDEL_SRC := $(JORDEL_SRC)                    		\
	$(DIR)/mod_gauss.f90                    	\
	$(DIR)/mod_basis_at_gaussian.f90        	\
	$(DIR)/hermite_1d.f90                   	\
	$(DIR)/initialise_basis.f90             	\
	$(DIR)/basis_functions.f90              	\
	$(DIR)/basis_functions1.f90             	\
	$(DIR)/basis_functions2.f90             	\
	$(DIR)/interp.f90                       	\
	$(DIR)/interp_RZ.f90

JORPOL_SRC := $(JORPOL_SRC)                    		\
	$(DIR)/mod_gauss.f90                    	\
	$(DIR)/mod_basis_at_gaussian.f90        	\
	$(DIR)/hermite_1d.f90                   	\
	$(DIR)/initialise_basis.f90             	\
	$(DIR)/basis_functions.f90              	\
	$(DIR)/basis_functions1.f90             	\
	$(DIR)/basis_functions2.f90             	\
	$(DIR)/interp.f90                       	\
	$(DIR)/interp_RZ.f90

ENBIGGEN_SRC := $(ENBIGGEN_SRC)               		\
	$(DIR)/mod_gauss.f90                    	\
	$(DIR)/mod_basis_at_gaussian.f90        	\
	$(DIR)/hermite_1d.f90                   	\
	$(DIR)/initialise_basis.f90             	\
	$(DIR)/basis_functions.f90              	\
	$(DIR)/basis_functions1.f90             	\
	$(DIR)/basis_functions2.f90             	\
	$(DIR)/interp.f90                       	\
	$(DIR)/interp_RZ.f90

JOREK2_TARGET2VTK_SRC := $(JOREK2_TARGET2VTK_SRC) 	\
        $(DIR)/mod_gauss.f90                      	\
        $(DIR)/mod_basis_at_gaussian.f90          	\
        $(DIR)/initialise_basis.f90               	\
        $(DIR)/basis_functions.f90                	\
        $(DIR)/basis_functions1.f90               	\
        $(DIR)/basis_functions2.f90               	\
        $(DIR)/interp.f90                         	\
        $(DIR)/interp_RZ.f90

JOREK2_POWERS_SRC := $(JOREK2_POWERS_SRC) 		\
        $(DIR)/mod_gauss.f90                      	\
        $(DIR)/mod_basis_at_gaussian.f90          	\
        $(DIR)/initialise_basis.f90               	\
        $(DIR)/basis_functions.f90                	\
        $(DIR)/basis_functions1.f90               	\
        $(DIR)/basis_functions2.f90               	\
        $(DIR)/interp.f90                         	\
        $(DIR)/interp_RZ.f90                      	\
        $(DIR)/hermite_1d.f90

JOREK2_IMPORT_PERTURBATION_SRC := $(JOREK2_IMPORT_PERTURBATION_SRC)	\
	$(DIR)/mod_gauss.f90                    	\
	$(DIR)/mod_basis_at_gaussian.f90        	\
	$(DIR)/initialise_basis.f90             	\
	$(DIR)/basis_functions2.f90 			\
	$(DIR)/basis_functions1.f90 			\
	$(DIR)/basis_functions.f90 			\
        $(DIR)/interp.f90                         

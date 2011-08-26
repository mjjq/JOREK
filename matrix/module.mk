DIR = matrix

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 			\
	$(DIR)/coicsr.f90 				\
	$(DIR)/construct_matrix.f90 			\
	$(DIR)/construct_matrix_murge.f90 		\
	$(DIR)/locate_irn_jcn.f90 			\
	$(DIR)/reduce.f90 				\
	$(DIR)/scale_global_matrix.f90 			\
	$(DIR)/sort_I_mrgrnk.f90 			\
	$(DIR)/global_matrix_structure.f90		\
	$(DIR)/mod_global_distributed_matrix.f90	\
	$(DIR)/ch_nod_rhs_elm.f90			\
	$(DIR)/ch_node_struct.f90
	
JOREK2_FOUR_SRC := $(JOREK2_FOUR_SRC) 			\
	$(DIR)/mod_global_distributed_matrix.f90	\

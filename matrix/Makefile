include ../config.in

all:  modules 			\
      coicsr.o 			\
      ch_nod_rhs_elm.o          \
      construct_matrix.o 	\
      construct_matrix_murge.o 	\
      locate_irn_jcn.o 		\
      reduce.o 			\
      scale_global_matrix.o 	\
      sort_I_mrgrnk.o 		\
      global_matrix_structure.o \
      ch_node_struct.o 

modules : mod_global_distributed_matrix.o

clean :
	rm -f *.o *.mod

include ../../config.in	

all: modules 				\
     initialise_parameters.o 		\
     sources.o 				\
     element_matrix.o 			\
     element_matrix_fft.o 		\
     apply_fft.o 			\
     boundary_conditions.o 		\
     initial_conditions.o		\
     boundary_matrix.o   		\
     boundary_matrix_open.o   		

modules : mod_parameters.o mod_phys_module.o mod_equation_variables.o 

clean :
	rm -f *.mod *.o

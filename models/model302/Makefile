include ../../config.in	

all: modules 				\
     initialise_parameters.o 		\
     sources.o 				\
     element_matrix.o 			\
     element_matrix_fft.o 		\
     boundary_matrix.o   		\
     boundary_matrix_open.o   		\
     boundary_conditions.o 		\
     boundary_conditions_murge.o 	\
     initial_conditions.o

modules : mod_parameters.o mod_phys_module.o 

clean :
	rm -f *.mod *.o

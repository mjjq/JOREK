include ../../config.in	

all: modules 				\
     initialise_parameters.o 		\
     element_matrix.o 			\
     boundary_conditions.o 		\
     initial_conditions.o

modules : mod_parameters.o mod_phys_module.o 

clean :
	rm -f *.mod *.o

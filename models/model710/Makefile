include ../../config.in	

all: modules 				\
     initialise_parameters.o 		\
     element_matrix.o 			\
     element_matrix_fft.o               \
     sources.o                          \
     boundary_conditions.o 		\
     initial_conditions.o

modules : mod_parameters.o mod_phys_module.o 

clean :
	rm -f *.mod *.o

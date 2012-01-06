include ../config.in	

all: modules \
     current.o \
     density.o \
     element_matrix_GS.o \
     element_matrix_GS_perturbation.o \
     element_matrix_GS_inverse.o \
     element_matrix_poisson.o \
     element_matrix_poisson_inverse.o \
     poisson.o \
     equilibrium.o \
     ffprime.o \
     chgmt_node.o \
     temperature.o \
     temperature_i.o \
     temperature_e.o \
     pellet_source.o \
     add_pellet.o    \
     log_parameters.o \
     read_num_profiles.o \
     derive_num_profiles.o \
     preset_parameters.o \
     initialise_and_broadcast_parameters.o

modules : mod_parameters.o mod_phys_module.o mod_pellet.o mod_nodes_elements.o mod_avoid_neg_dens.o mod_domains.o

clean :
	rm -f *.mod *.o

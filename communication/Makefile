
include ../config.in

all: mod_phys_module.o \
     mod_mumps.o \
     broadcast_elements.o \
     broadcast_boundary.o \
     broadcast_nodes.o \
     broadcast_phys.o \
     broadcast_num_profiles.o \
     distribute_harmonics.o \
     distribute_nodes_elements.o \
     distribute_vector.o  \
     export_grid.o \
     export_pov.o \
     export_restart.o \
     import_restart.o \
     update_deltas.o \
     update_values.o \
     export_boundary.o \
     vertex_is_local.o \
     export_helena.o  \
     split_broadcast.o

modules : mod_live_data.o

clean :
	rm -f *.o *.mod

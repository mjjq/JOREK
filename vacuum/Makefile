include ../config.in	

all : modules                    \
	ideal_wall.o                 \
	ideal_wall_starwall.o        \
	resistive_wall_starwall.o    \
	vacuum_old.o
      
modules : mod_vacuum.o           \
	mod_vacuum_response.o        \
	mod_vacuum_equilibrium.o
    
clean :
	rm -f *.o *.mod

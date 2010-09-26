include ../config.in	

all :   modules                   \
	get_vacuum_response.o     \
	ideal_wall.o              \
	ideal_wall_starwall.o     \
	resistive_wall_starwall.o \
	vacuum.o                  \
	vacuum_equil.o
      
modules : mod_vacuum_response.o     
    
clean :
	rm -f *.o *.mod
    

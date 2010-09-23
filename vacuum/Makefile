include ../config.in	

all :   modules       \
	ideal_wall.o      \
	resistive_wall.o  \
	vacuum.o          \
	vacuum_equil.o
      
modules : mod_vacuum_response.o     
    
clean :
	rm -f  *.o
    

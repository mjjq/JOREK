include ../config.in

.SUFFIXES: .c
%.o:%.c
	$(CC) $(CFLAGS)           -c $< -o $@  $(INCLUDES)

all : fortran_pthread.o 

modules : 

clean :
	rm -f *.mod *.o

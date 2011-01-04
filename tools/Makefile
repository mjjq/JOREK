include ../config.in

.SUFFIXES: .c
%.o:%.c
	$(CC) $(CFLAGS)           -c $< -o $@  $(INCLUDES)

all : fortran_pthread.o

fortran_pthread.c

include ../config.in	

all: modules \
     basis_functions2.o \
     basis_functions1.o \
     basis_functions.o \
     bezier_1d.o \
     hermite_1d.o \
     hermite_elements.o \
     initialise_basis.o \
     interp.o \
     interp_RZ.o
     
     
modules : mod_gauss.o \
          mod_basis_at_gaussian.o 

clean :
	rm *.mod *.o

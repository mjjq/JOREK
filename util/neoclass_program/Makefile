FC =ifort
#FC ?=gfortran
FFLAGS?=-g -C

#INCLUDES = -I$(MKL_HOME)/include

BLAS_HOME = ${MKLROOT}/lib/intel64
BLASLIB =  $(BLAS_HOME)/libmkl_intel_lp64.a $(BLAS_HOME)/libmkl_sequential.a $(BLAS_HOME)/libmkl_core.a 
LFLAGS =$(FFLAGS) #-L$(MKL_HOME)/lib
LIBS = $(BLASLIB)

FOBJS = countlines.o calculs_coefs.o neoclassical_program.o 

neoclassical_program : $(FOBJS)
	$(FC) $(INCLUDES) -o $@ $^ $(LFLAGS) $(LIBS)

%.o : %.f90
	$(FC) $(INCLUDES)-c $(FFLAGS) $<

#countlines.o : countlines.f90
#	$(FC) -c $(FFLAGS) -o $@ $<

neoclassical_program.o : calculs_coefs.o countlines.o

clean : 
	rm -f *.o *.mod

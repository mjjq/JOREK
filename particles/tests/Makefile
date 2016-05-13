INCLUDES=-I. -I.. ../openadas.o
INCLUDES_MC=$(INCLUDES) ../mod_ionisation_recombination.o ../../elements/interp_PRZ.o ../../models/mod_phys_module.o ../../elements/basis_functions3.o
FFLAGS=-fbounds-check -Wall -llapack -g

all: test_mc_time test_coronal_time test_l2d_interp test_linear_spline

test_mc_time: test_mc_time.f90 Makefile ../openadas.o ../mod_ionisation_recombination.o
	      gfortran $(INCLUDES_MC) $< $(FFLAGS) -o $@

test_coronal_time: test_coronal_time.f90 Makefile ../openadas.o ../mod_coronal.o
	gfortran $(INCLUDES) ../mod_coronal.o $< $(FFLAGS) -o $@

test_l2d_interp: test_l2d_interp.f90 Makefile ../openadas.o
	gfortran -fcheck=all -ffpe-trap=zero,overflow,underflow $(INCLUDES) $< $(FFLAGS) -o $@

test_linear_spline: test_linear_spline.f90 Makefile ../openadas.o
	gfortran -fcheck=all -ffpe-trap=zero,overflow,underflow $(INCLUDES) $< $(FFLAGS) -o $@

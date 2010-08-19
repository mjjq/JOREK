

plots/ppplib.o: plots/ppplib.f
	$(FC) $(FFLAGS_NOBOUNDS)  -c $< -o $@  $(INCLUDES)

solvers/dPackgmres.o: solvers/dPackgmres.f
	$(FC) $(FFLAGS_FIXEDFORM) -c $< -o $@  $(INCLUDES)

%.o: %.f90
	$(FC) $(FFLAGS)           -c $< -o $@  $(INCLUDES)

%.o: %.f
	$(FC) $(FFLAGS)           -c $< -o $@  $(INCLUDES)


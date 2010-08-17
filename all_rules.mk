

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(INCLUDES)

%.o: %.f
	$(FC) -c $< -o $@ $(FFLAGS) $(INCLUDES)
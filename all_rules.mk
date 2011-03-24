

plots/ppplib.o: plots/ppplib.f
	@echo "*** $< ***"
	$(FC) $(FFLAGS_NOBOUNDS) -c $< -o $@  $(INCLUDES)

solvers/dPackgmres.o: solvers/dPackgmres.f
	@echo "*** $< ***"
	$(FC) $(FFLAGS_FIXEDFORM) -c $< -o $@  $(INCLUDES)

%.o: %.f90
	@echo "*** $< ***"
	$(FC) $(FFLAGS) -c $< -o $@  $(INCLUDES)

%.o: %.f
	@echo "*** $< ***"
	$(FC) $(FFLAGS) -c $< -o $@  $(INCLUDES)

%.o: %.c
	@echo "*** $< ***"
	$(CC) $(CFLAGS) -c $< -o $@  $(INCLUDES)

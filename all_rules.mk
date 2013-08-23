plots/ppplib.o: plots/ppplib.f Makefile Makefile.inc
	@echo "*** $< ***"
	$(FC) $(FFLAGS_NOBOUNDS) -c $< -o $@  $(INCLUDES)

solvers/dPackgmres.o: solvers/dPackgmres.f Makefile Makefile.inc
	@echo "*** $< ***"
	$(FC) $(FFLAGS_FIXEDFORM) -c $< -o $@  $(INCLUDES)

%.o: %.f90 Makefile Makefile.inc
	@echo "*** $< ***"
	$(FC) $(FFLAGS) -c $< -o $@  $(INCLUDES)

%.o: %.f Makefile Makefile.inc
	@echo "*** $< ***"
	$(FC) $(FFLAGS) -c $< -o $@  $(INCLUDES)

%.o: %.c Makefile Makefile.inc
	@echo "*** $< ***"
	$(CC) $(CFLAGS) -c $< -o $@  $(INCLUDES)

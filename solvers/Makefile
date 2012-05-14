include ../config.in	

all: 	modules 		\
	gmres_driver.o 		\
	gmres_matrix_vector.o 	\
	gmres_precondition.o 	\
	dPackgmres.o  		\
	initialise_mumps.o 	\
	initialise_pastix.o 	\
	mnewtax.o 		\
	root.o 			\
	solve_M2.o 		\
	solve_matrix_n.o 	\
	solve_mumps_all.o 	\
	solve_pastix_all.o 	\
	solvP3.o 		\
	solve_murge_all.o

modules : 	mod_mumps.o 	\
		mod_pastix.o	\
		mod_murge.o
clean :
	rm -f *.mod *.o	  

# Time evolution of coronal equilibrium
In the coronal equilibrium calculation, the equilibrium is reached quite quickly for W.
Setup the situation with a uniform distribution over each charge state.
Using backward euler, the system will equilibriate in one timestep if this step is > roughly 10^-6.
This can be tested with the program test_coronal_time.f90
(compile with gfortran test_coronal_time.f90 ../openadas.o ../mod_coronal.o -llapack -o test_coronal_time after compiling jorek2_particles)

See images for more results. Noise on the low-q side of the graph is probably caused by precision of the coefficients at these high temperatures / numerical issues as they show up after some time only.

In the future we will therefore perform one single timestep with the backward euler method (theta = 1) to find the charge state


# Evaluation of the coronal equilibrium
The program test_coronal_eq can be used to calculate the coronal equilibrium for several densities and temperatures.


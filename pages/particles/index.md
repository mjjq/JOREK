title: JOREK2_particles

## JOREK2_particles is a MC kinetic/guiding-center particle transport code

## Features
* Multiple particle pushers available
    - Boris method (kinetic)
    - Guiding-center
    - More planned
* [Integration with JOREK](jorek_integration.html)
    - Feedforward integration for EM fields during MHD (in)stabilities
    - Feedback to JOREK planned
* Written in modern fortran 2003
* Particle initialization with PCG random numbers or Sobol sequences.
* Projection of particle positions onto JOREK finite elements
* Parallelization with MPI and OpenMP
* Load-balancing
* Parallel HDF5 input/output
* Diagnostics framework
* Ionization/recombination based on the ADAS ADF11 database
* Impurity radiation diagnostic (in progress)
* Particle-background collisions based on a binary collision method (in progress)


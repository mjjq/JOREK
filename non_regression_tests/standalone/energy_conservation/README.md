# Particle conservation laws test
This testcase tracks particles in a grad-shafranov equilibrium as calculated by JOREK and checks conservation of
- energy
- generalized toroidal momentum
of particles released at random positions in the domain.

## Conservation of energy
The energy of the particle is given by the kinetic energy and the electric potential energy.
This is calculated from the potential $U$ and the $F$-function.
\[ E = 0.5 m_p \mathbf{v}^2 + 0.5 q / m_p * F0 * U \]

In the JOREK formulation the positions and fields are known at the full timesteps and the velocities at the half-timesteps.
We must also convert the equation for the energy to JOREK units.
This leads to the following expression to calculate the energy at timestep $n$, (where $v = 0.5 (v^n+v^{n+1})$ and the dot product distributes).
``
E = 0.25 * particle%mass * ATOMIC_MASS_UNIT * ( dot_product(v,v) + dot_product(v_old, v_old)) 
  + particle%q * EL_CHG * t_norm * F0 * U
``

## Conservation of momentum
In this toroidal field we are actually concerned with the conservation of canonical momentum.
This boils down to the addition of a term with the local magnetic flux.
\[ p = r v_\phi - q/m \psi_{} \]
In JOREK the fields are not normalized so we need not alter this equation.
It is still necessary to obtain the time-centered velocities instead of the velocities at the half-steps.
``
mp = x(1) * v(3) + qom * psi
``
where `qom` is the particle charge over particle mass.


## Testcase setup
Run the testcase with `run.sh` to locally run a particle job.
See `run.sh -h` for help with the command.

### Running on a cluster
Submit the testcase job as follows: (compile jorek2_particles beforehand) (omp only supported, so at most 1 node)
qsub ./run.sh -q ib_gen8 -l nodes=1:ppn=16 -v file=jet_equil_short
Run all jobs with:
for file in jet*; do qsub ./run.sh -q ib_gen8 -l nodes=1:ppn=16 -v file=$file; done

### Dependencies
Install numpy and matplotlib on your system (anaconda is the simplest way)

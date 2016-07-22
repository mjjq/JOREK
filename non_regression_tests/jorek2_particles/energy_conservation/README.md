# Particle conservation laws test
This testcase tracks particles in a field as calculated by JOREK and checks conservation of
- energy
- generalized toroidal momentum
of particles released at random positions in the domain.

## Conservation of energy
The energy of the particle is given by the kinetic energy and the electric potential energy.
This is calculated from the potential $U$ and the $F$-function.
\[ E = 0.5 m_p \mathbf{v}^2 + 0.5 q / m_p * F0 * U \]

In the JOREK formulation the positions and fields are known at the full timesteps and the velocities at the half-timesteps.
We must also convert the equation for the energy to JOREK units.
This leads to the following expression to calculate the energy at timestep $n$.
``
E = 0.5 * particle%mass * ATOMIC_MASS_UNIT * dot_product(v,v)
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


## Results
Output will be many energies and momenta at all nout timesteps, plots of the distribution, traces of the output.
Combine many of these (for a timestep scan) with the dt_analyze.sh script, which produces a graph of the error vs the timestep.
It is important that the total simulation time is constant, because some effects scale as the number of steps and other with the square of the timestep size.

### Energy conservation
The Boris method integrator should conserve energy exactly. Looking at different integration lengths we see that this is true down to the floating point error, which increases with the square root of the number of steps.
The normalized traces of particle energies over many timesteps show no drift. Since these are normalized the change induced by an error of order $\epsilon$ varies with the initial energy, producing the plot below.
![Energy versus step number](jet_eq_t0.00015_out/energy_trace.png)

We can look at the growth rate of the error in energy, which is shown below. The variance should scale with the square root of the number of steps.
![Growth rate of energy. The straight line is a square-root of the number of timesteps and serves as a guide to the eye.](energy_growth.png)
![Energy variance. For a random walk this scales with the square root of the number of steps, which corresponds nicely to the line plotted in the figure.](energy_variance.png)

### Momentum conservation
The toroidal momentum consists of a velocity-component and an position-component, where the error in the position-component dominates.
If we look at a time-trace, it can be seen that there is a lot of variation, but because the velocity is reproduced well in the Boris scheme this is counteracted on the two halves of the orbit.
![Normalized momentum versus step number](jet_eq_t0.00015_out/momentum_trace.png)

This is also exemplified by the scaling with the timestep size, which corresponds to the error scaling of the Boris method (2nd order).
![Growth rate of momentum.](momentum_growth.png)
The variance is also reduced by a second-order scaling with the timestep.
![Momentum variance](momentum_variance.png)

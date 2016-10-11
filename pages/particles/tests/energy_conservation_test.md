title: Energy conservation test

This testcase tracks particles in a field as calculated by JOREK and
checks conservation of energy and generalized toroidal momentum of
particles released at random positions in the domain.

## Conservation of energy

The energy of the particle is given by the kinetic energy and the
electric potential energy. This is calculated from the potential \(U\)
and the \(F\)-function.
\[ E = 0.5 m_p^2 + 0.5 q / m_p * F_0 * U \]

In the JOREK formulation the positions and fields are known at the full
timesteps and the velocities at the half-timesteps. We must also convert
the equation for the energy to JOREK units. This leads to the following
expression to calculate the energy at timestep \(n\).
```fortran
E = 0.5 * particle%mass * ATOMIC_MASS_UNIT * dot_product(v,v) + particle%q * EL_CHG * t_norm * F0 * U
```

## Conservation of momentum

In this toroidal field we are also concerned with the conservation
of canonical momentum. This boils down to the addition of a term with
the local magnetic flux. \( p = r v_\phi - q/m \).
In JOREK the fields are not normalized so we need not alter this equation. It is
still necessary to obtain the time-centered velocities instead of the
velocities at the half-steps. 
```fortran
mp = x(1) * v(3) + qom * psi
```
where qom is the particle charge over particle mass.

## Testcase setup

Run the testcase with `sh run.sh` to locally run a particle job.

### Running on a cluster

Submit the testcase job as follows: (compile jorek2\_particles
beforehand) (only OMP supported, so at most 1 node)
```bash
qsub ./run.sh -q ib_gen8 -l nodes=1:ppn=16 -v file=jet_equil_short
```
Run all jobs with:
```bash
for file in jet*; do
  qsub ./run.sh -q ib_gen8 -l nodes=1:ppn=16 -v file=$file;
done
```

### Dependencies

Install numpy and matplotlib on your system (anaconda is the simplest
way)

### Results

Output will be many energies and momenta at all nout timesteps, plots of
the distribution, traces of the output. Combine many of these (for a
timestep scan) with the dt\_analyze.sh script, which produces a graph of
the error vs the timestep. It is important that the total simulation
time is constant, because some effects scale as the number of steps and
other with the square of the timestep size. In these cases we have
simulated for a total time of 150.0 JOREK units, with timesteps of
0.00015 to 0.060. The number of steps is thus at most \(10^6\).

#### Energy conservation

The Boris method integrator should conserve energy exactly. Looking at
different integration lengths we see that this is true down to the
floating point error, which increases with the square root of the number
of steps. The normalized traces of particle energies over many timesteps
show no drift. Since these are normalized the change induced by a
numerical error of size \(\epsilon\) varies with the initial energy,
producing the plot below.

![energy traces](energy_conservation_jet_eq_t0.00015_out_energy_trace.png)

We can look at the growth rate of the error in energy, which is shown
below. The variance should scale with the square root of the number of
steps.

![energy growth rate](energy_conservation_energy_growth.png)![Energy variance](energy_conservation_energy_variance.png)

#### Momentum conservation

The toroidal momentum consists of a velocity-component and an
position-component, where the error in the position-component dominates.
If we look at a time-trace, it can be seen that there is a lot of
variation, but because the velocity is reproduced well in the Boris
scheme this is counteracted on the two halves of the orbit.

![energy traces](energy_conservation_jet_eq_t0.00015_out_momentum_trace.png)

This is also exemplified by the scaling with the timestep size, which
corresponds to the error scaling of the Boris method (2nd order).
The variance is also reduced by a second-order scaling with the
timestep.

![energy growth rate](energy_conservation_momentum_growth.png)![Energy variance](energy_conservation_momentum_variance.png)

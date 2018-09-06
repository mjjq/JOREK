
# The Penning trap

The Penning trap consists of a constant, uniform magnetic field
\(\mathbf{B}_0\) in the \(z\)-direction, and a quadrupole electric
field. The configuration is shown below. (Image from [Wikimedia, by Arian
Kriesch](https://commons.wikimedia.org/wiki/File:Penning_Trap.svg))

![penning trap](|media|/tests/penning/500px-Penning_Trap.png)


## Parameters

The independent parameters used in the test are

| Parameter | Value | Description |
|-----------|-------|-------------|
| \(\omega_e\) | 4.9 rad/s | Oscillation frequency of the particle due to the electric field |
| \(\omega_b\) | 25.0 rad/s | Oscillation frequency of the particle due to the magnetic field |
| \(\epsilon\) | -1 | Polarisation of the penning trap | 
| \(\mathbf{x}_0\) | (10,0,0) m | Initial position (in xyz coordinates) |
| \(\mathbf{v}_0\) | (50,0,20) m/s | Initial velocity (in xyz coordinates) |

For these parameters the particle trajectory has been [calculated](|media|/tests/penning/penning_trap.nb).
The x-y position of the particle over time is shown in the figure below.
The motion in the \(z\)-direction is decoupled and is a simple harmonic oscillator.

![xy-trajectory](|media|/tests/penning/penning_xy.png)![z-trajectory](|media|/tests/penning/penning_z.png)

## Results
The trajectory of a particle is followed from \(\mathbf{x}_0\) at \(t=0 s\) to \(t=16 s\).
An example trajectory, calculated with the [[mod_boris]] method and relatively large timesteps (\(\delta t = 0.01\)) is shown below.

![xy-trajectory-boris](|media|/tests/penning/penning_xy_boris.png)

### Comparing pushers
The different pushers in the [[pusher_test]] perform as follows on this case

![pusher-test-penning](|media|/tests/all_pushers/penning.png)
![pusher-test-penning](|media|/tests/all_pushers/penning_time.png)

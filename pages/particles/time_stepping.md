title: Time stepping

# Time stepping in J2P
The timestep of each particle group is calculated before the simulation based on
the requested timestep and the requested sync times.

Sync times are calculated based on each event, which can request sync of any or
all of the groups' times. Imagine a simulation export or diagnostic tool that
requires the values at that time.

Each of the event intervals can be specified as a $\Delta T_i$ for event $i$.
An optional weight can be set, or the default for this type of event will be
used. The actual timesteps and event intervals $\underline x$ will be calculated to minimize
a cost function
\(
  ||c - Ax||^2
\)



The base timestep, or tick, used in the simulation is denoted
\(
  \Delta t.
\)
All other times used in the simulation must be integer multiples of this tick.
There are usually at least two timesteps imposed upon the simulation.
These are the input file timestep and the steps from any outputs.
Let's call these $\Delta T_i$.
Together, these determine the fundamental time unit of the simulation, as the
greatest common divisor.
This is however only defined for integers. We can calculate an approximation if
it is allowed that some timesteps change slightly. 

Each of the $\Delta T_i$ has an associated allowed relative change $\epsilon_i$
which is $1\cdot10^{-5}$ by default, but can be set by the user.


\(
  \Delta T = 
\)









# old
The time-scale is determined by the base timestep, set in the input parameters.
Let's call this variable
\(
  \Delta T
\)

All other events then need to be expressed in integer multiples or divisions of this time.
We will use the 

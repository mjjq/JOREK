subroutine pellet_source(pellet_amplitude,pellet_R,pellet_Z,pellet_phi,pellet_radius, &
                         pellet_sig,pellet_length,R,Z,phi,particle_source)

implicit none

real*8 :: R, Z, phi, particle_source
real*8 :: pellet_amplitude, pellet_R, pellet_Z, pellet_phi, pellet_radius, pellet_sig, pellet_length
real*8 :: radius, atn, atn_phi, PI

PI = 3.14159265358979

if (phi .gt. PI) phi = 2*PI - phi

radius = sqrt((R-pellet_R)**2 + (Z-pellet_Z)**2)

atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))

atn_phi = (0.5d0 - 0.5d0*tanh(abs((phi- pellet_phi))/pellet_length))

particle_source = pellet_amplitude * atn * atn_phi

return
end

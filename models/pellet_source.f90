subroutine pellet_source(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                         pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
			 R,Z,psi,phi,particle_source)

use constants

implicit none

real*8 :: R, Z, psi, phi, particle_source
real*8 :: pellet_amplitude, pellet_R, pellet_Z, pellet_phi, pellet_radius, pellet_sig, pellet_length
real*8 :: pellet_psi, pellet_delta_psi, radius, atn, atn_psi, atn_phi

if (phi .gt. PI) phi = 2*PI - phi

radius = sqrt((R-pellet_R)**2 + (Z-pellet_Z)**2)

atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))

atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))

atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))

particle_source = pellet_amplitude * atn * atn_phi * atn_psi

return
end

module mgi_module

use constants

real*8 :: total_n_particles_inj
real*8 :: total_n_particles_plasma
real*8 :: total_n_particles_inj_all

contains 

subroutine mgi_source(mgi_amplitude,mgi_R,mgi_Z,mgi_phi,mgi_radius, &
                         mgi_sig,mgi_length,R,Z,phi,particle_source)

implicit none

real*8 :: R, Z, phi, particle_source
real*8 :: mgi_amplitude, mgi_R, mgi_Z, mgi_phi, mgi_radius, mgi_sig, mgi_length
real*8 :: radius, atn, atn_phi, PI

PI = 3.14159265358979

if (phi .gt. PI) phi = 2*PI - phi

radius = sqrt((R-mgi_R)**2 + (Z-mgi_Z)**2)

atn     = (0.5d0 - 0.5d0*tanh((radius - mgi_radius)/mgi_sig))

atn_phi = (0.5d0 - 0.5d0*tanh((phi- mgi_phi)/mgi_length))

particle_source = mgi_amplitude * atn * atn_phi

return
end subroutine mgi_source

subroutine update_mgi(my_id,node_list,element_list)

use constants
use data_structure
use phys_module
use mpi_mod

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer :: my_id, ierr
real*8  :: density, density_in, density_out, pressure, pressure_in,pressure_out

call Integrals_3D(my_id, node_list, element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

total_n_particles_inj = total_n_particles_inj * n_zero * tstep
total_n_particles_inj_all = total_n_particles_inj_all + total_n_particles_inj
total_n_particles_plasma = total_n_particles_plasma * n_zero

if (my_id .eq. 0) then
  write(*,*) 'total neutrals particles injected at that timestep =', total_n_particles_inj
  write(*,*) 'total neutrals particles into the plasma =', total_n_particles_plasma
  write(*,*) 'total neutrals particles injected since the start of the simulation = ', total_n_particles_inj_all
  write(*,*) 'Check of density conservation'
endif

end subroutine update_mgi

end module mgi_module

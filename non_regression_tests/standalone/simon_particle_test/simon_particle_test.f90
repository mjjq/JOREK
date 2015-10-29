!> This program tests the reproduction of particle trajectories in a
!! large aspect ratio tokamak with circular flux surfaces.
!!
!! See http://www2.ipp.mpg.de/~Simon.Pinches/thesis/node57.html for details.
!! This program reads a namelist input file containing the JOREK settings, and uses
!! tstep and nstep to determine the integration parameters.
program simon_particle_test

use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements
use mod_particles
use clock_module
use parameters
use constants

implicit none

! Hard-coded parameters
integer, parameter :: my_id=0

! Internal variables
integer :: ierr,i,j

type (type_particle_list) :: particle_list

write(*,*) '***************************************'
write(*,*) '* JOREK2 : Simon Particle test        *'
write(*,*) '***************************************'

!! Read input parameters
call initialise_parameters(my_id, "__NO_FILENAME__")

!! Read restart file
call import_binary_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)

!! Setup parameters and grid
call initialise_basis                              ! define the basis functions at the Gaussian points
call update_neighbours(element_list,node_list)     ! update neighbour information in the element_list

!! Initialize particles
call initialise_particles_simon(node_list,element_list,particle_list)

!! Perform time-stepping manually
do i=1,nstep
  call update_particles(my_id, particle_list, tstep, 1, 0.d0)
  !! Save particle position to file
  do j=1,3
    write(*,*) j,particle_list%particle(j)%x
  enddo
enddo


contains

!> Setup the 3 particles we will use and test some invariants
subroutine initialise_particles_simon(node_list,element_list,particle_list)

use mod_particles
use constants
use data_structure
use phys_module, only: F0, central_density, central_mass

implicit none

type (type_node_list), intent(in)        :: node_list
type (type_element_list), intent(in)     :: element_list
type (type_particle_list), intent(inout) :: particle_list

type (type_particle)      :: particle

real*8  :: particle_energy(3), particle_energy_perp(3), R_in, Z_in, phi_in, R_out, Z_out, v_R, v_phi
real*8  :: P(1), P_s(1), P_t(1), R, R_s, R_t, Z, Z_s, Z_t, B_field(3), B_0, st_jac, v_norm, v_perp2, v_par
real*8  :: psi_s, psi_t, psi_R, psi_Z, s_elm, t_elm, mass_ion
integer :: i_var(1), i_elm, ifail, i_part


particle_list%n_particles = 3
allocate(particle_list%particle(particle_list%n_particles))

! Start position of all 3 particles
R_in   = 3.025
Z_in   = 0.0
phi_in = 0.0

call find_RZ(node_list,element_list,R_in,Z_in,R_out,Z_out,i_elm,s_elm,t_elm,ifail)

particle_energy      = (/ 170., 50.,    164.    /)      ! [eV]
particle_energy_perp = (/  40., 49.585, 161.832 /)      ! [eV]

mass_ion = mass_proton * central_mass
v_norm = sqrt(mu_zero * mass_ion * central_density * 1.d20)

i_var = (/ 1 /)
call interp_PRZ(node_list,element_list,i_elm,i_var,1,s_elm,t_elm,phi_in,P,P_s,P_t,R,R_s,R_t,Z,Z_s,Z_t)

st_jac = R_s * Z_t - R_t * Z_s
psi_s  = P_s(1); psi_t = P_t(1);
psi_R  = (  psi_s * Z_t - psi_t * Z_s ) / st_jac
psi_Z  = (- psi_s * R_t + psi_t * R_s ) / st_jac

B_field     = (/ + psi_Z, - psi_R, F0 /) / R
B_0         = sqrt(dot_product(B_field,B_field))

do i_part=1,3
  v_par = sqrt((particle_energy(i_part) - particle_energy_perp(i_part))* 2. / mass_ion * el_chg)
  v_phi = v_par * B_0 / B_field(3)
  v_perp2 = particle_energy_perp(i_part)* 2. / mass_ion * el_chg

  write(*,'(A,3e16.8)') ' perpendicular velocity : ',sqrt(v_perp2),v_norm,sqrt(v_perp2) * v_norm
  write(*,'(A,3e16.8)') ' gyro radius            : ',mass_ion * sqrt(v_perp2) / (el_chg * B_0)
  write(*,'(A,3e16.8)') ' gyro frequency         : ',el_chg * B_0 / mass_ion, el_chg * B_0 / mass_ion * v_norm

  v_R = sqrt(v_perp2 - v_phi**2*(1. - (B_field(3)/B_0)**2)**2 - v_phi * B_field(3)**2 * B_field(2)**2 / B_0**4)

  write(*,'(A,3e16.8)') 'CHECK energy : ', v_R,v_phi,(v_R**2+v_phi**2) * 0.5 * mass_ion / el_chg

  particle%v(1) = - v_R * v_norm
  particle%v(2) = 0.
  particle%v(3) = - v_phi * v_norm

  particle%x       = (/ R_out, Z_out, 0. /)        !< particle position in real space
  particle%st      = (/ s_elm, t_elm /)            !< particle position in the finite element (i_elm)
  particle%i_elm   = i_elm
  particle%q       = 1                             !< charge
  particle%mass    = central_mass                  !< mass
  particle%weight  = 1.                            !< weight (i.e. number of particles)

  particle_list%particle(i_part) = particle

  write(*,'(A,3e16.8)') 'Velocity : ',particle%v
enddo
end subroutine initialise_particles_simon
end program simon_particle_test

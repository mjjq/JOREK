module pellet_module


real*8 :: total_pellet_particles   !< the (total) pellet particles added in this timestep
real*8 :: total_plasma_particles   !< the total plasma density (before this timestep)
real*8 :: total_pellet_volume      !< the volume of the simulated pellet in this timestep

real*8 :: phys_pellet_volume       !< the physical pellet radius (in m^3)
real*8 :: pellet_volume            !< approximated value of simulated pellet volume
real*8 :: pellet_atomic            !< atomic number of pellet mass

contains

subroutine pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                          pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
                          R,Z,psi,phi, r0, T0, central_density, pellet_particles, pellet_density, &
                          particle_source, volume_source)

implicit none

!input variables
real*8 :: R, Z, psi, phi            ! position whereh the particle source will be calculated
real*8 :: T0, r0                    ! local temperature and mass density (JOREK normalised)
real*8 :: central_density           !< central plasma density (in units 10^20 m^-3)
real*8 :: pellet_particles          !< total number of particles in the pellet
real*8 :: pellet_density            !< pellet density (units 10^20 m^-3)
real*8 :: pellet_amplitude          !< amplitude of paricle source (when not using ablation model)
real*8 :: pellet_R, pellet_Z        !< position of the pellet (phi=0)
real*8 :: pellet_phi                !< length of pellet in toroidal direction
real*8 :: pellet_radius             !< pellet size (radius) in poloidal plane
real*8 :: pellet_sig, pellet_length !< sigmas of pellet source in poloidal and toroidal direction
real*8 :: pellet_psi, pellet_delta_psi

!output variables
real*8 :: particle_source           !< particle source (JOREK normalised units)
real*8 :: volume_source             !< volume of the pellet source (variable used to integrate total pellet volume)

!local variables
real*8  :: PI, radius, atn, atn_psi, atn_phi, atomic_mass, ablation_rate

particle_source = 0.d0
volume_source   = 0.d0

PI = 2.d0 *asin(1.d0)

if (pellet_amplitude .gt. 0.) then             ! use the fixed source pellet model 

  if (phi .gt. PI) phi = 2*PI - phi

  radius = sqrt((R-pellet_R)**2 + (Z-pellet_Z)**2)

  atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
  atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
  atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))

  particle_source = pellet_amplitude * atn * atn_phi * atn_psi

!S.F. modified here for introducing moving pellet...
else if (pellet_particles .gt. 0.) then
 
  if (phi .gt. PI) phi = 2*PI - phi

  radius = sqrt((R-pellet_R)**2 + (Z-pellet_Z)**2)
  atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
  atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))
  atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
  
  pellet_volume = PI * pellet_radius**2 * pellet_R * pellet_phi ! simulated pellet volume
    
  phys_pellet_volume = pellet_particles /pellet_density         ! physical pellet volume 

! the number of particles ablated from the physical pellet in units 10^20 m^-3 per unit of JOREK time

  pellet_atomic = 2.d0

  ablation_rate = 1.62d5 * central_density**(-0.77) * T0**(1.72) * r0**(0.45) * phys_pellet_volume**(0.48) * pellet_atomic**0.217
  
! particle source in JOREK normalisation

  particle_source = ablation_rate / central_density  * atn * atn_phi / pellet_volume
  
  volume_source   = atn * atn_phi

end if

return
end subroutine pellet_source2

subroutine update_pellet(my_id,node_List,element_list)
!******************************************************************************
! routine updates the pellet position and the size of the simulated           *
! and physical pellet sizes (from the integral of the pellet particle source) *
!******************************************************************************

use data_structure
use phys_module
implicit none

include 'mpif.h'

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer :: my_id, ierr
real*8  :: PI, zmu0, V_normalisation, density, density_in, density_out, pressure,pressure_in,pressure_out

PI   = 2.d0  * asin(1.d0)
zmu0 = 4.d-7 * PI

if (pellet_amplitude .gt. 0) return

! S.F. introduced pellet velocity here. 14/02/2011.

call Integrals_3D(my_id, node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

V_normalisation = 1.d0 / sqrt(central_density * 1.66d-7 * 2.d0 * zmu0) ! assumes Deuterium!

pellet_R = pellet_R + pellet_velocity_R * tstep / V_normalisation
pellet_Z = pellet_Z + pellet_velocity_Z * tstep / V_normalisation

PI = 2.d0*asin(1.d0)
 
total_pellet_particles = total_pellet_particles * central_density * tstep  ! check tstep

if (my_id .eq. 0) then

!  if (pellet_particles .gt. 0.0) then

    pellet_particles = max(pellet_particles - total_pellet_particles, 0.d0)
    
    write(*,'(i3,A,4e14.6)') my_id,' pellet (R,Z) =', pellet_R, pellet_Z,pellet_velocity_R/V_normalisation,pellet_velocity_Z/V_normalisation
    write(*,'(i3,A,4e14.6)') my_id,' total particles added in this step : ', pellet_R, total_pellet_particles, total_plasma_particles
    write(*,'(i3,A,4e14.6)') my_id,' remaining particles in pellet      : ', pellet_particles
    write(*,'(i3,A,4e14.6)') my_id,' pellet volume (sim,phys)           : ', total_pellet_volume,pellet_particles/pellet_density
   
!  end if

else 
  pellet_particles = 0.0
end if
 
call MPI_Bcast(pellet_particles,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

return 
 
end subroutine update_pellet

end module pellet_module

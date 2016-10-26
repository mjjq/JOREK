module pellet_module

use constants

real*8 :: total_pellet_particles   !< the (total) pellet particles added in this timestep
real*8 :: total_plasma_particles   !< the total plasma density (before this timestep)
real*8 :: total_pellet_volume      !< the volume of the simulated pellet in this timestep

real*8 :: phys_pellet_volume       !< the physical pellet radius (in m^3)
real*8 :: pellet_volume            !< approximated value of simulated pellet volume
real*8 :: pellet_atomic            !< atomic number of pellet mass

real*8 :: phys_ablation            !< physical ablation rate (non normalised)



real*8, allocatable  :: xtime_pellet_R(:)
real*8, allocatable  :: xtime_pellet_Z(:)
real*8, allocatable  :: xtime_pellet_psi(:)
real*8, allocatable  :: xtime_pellet_particles(:)
real*8, allocatable  :: xtime_phys_ablation(:)

contains

subroutine pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                          pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                          R,Z,psi,phi, r0, T0, central_density, pellet_particles, pellet_density, pellet_volume, &
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
real*8 :: pellet_ellipse            !< ellipticity of the pellet source in the poloidal plane
real*8 :: pellet_theta              !< orientation of the pellet ellipse in the poloidal plane
real*8 :: pellet_psi, pellet_delta_psi
real*8 :: pellet_volume

!output variables
real*8 :: particle_source           !< particle source (JOREK normalised units)
real*8 :: volume_source             !< volume of the pellet source (variable used to integrate total pellet volume)

!local variables
real*8  :: radius, atn, atn_psi, atn_phi, atomic_mass, ablation_rate

particle_source = 0.d0
volume_source   = 0.d0

if (pellet_amplitude .gt. 0.) then             ! use the fixed source pellet model 

pellet_particles = 0.0

  if (phi .gt. PI) phi = 2*PI - phi

  radius = sqrt(  (cos(pellet_theta)**2 + 1./pellet_ellipse**2 * sin(pellet_theta)**2)*(R-pellet_R)**2  &
                + (sin(pellet_theta)**2 + 1./pellet_ellipse**2 * cos(pellet_theta)**2)*(Z-pellet_Z)**2  &
                + 2.*(R-pellet_R)*(Z-pellet_Z)*sin(pellet_theta)*cos(pellet_theta) * (1./pellet_ellipse**2 - 1.) )

  atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
  atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
  atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))

  particle_source = pellet_amplitude * atn * atn_phi * atn_psi

!S.F. modified here for introducing moving pellet...
else if (pellet_particles .gt. 0.) then

  if (phi .gt. PI) phi = 2*PI - phi

  radius = sqrt(  (cos(pellet_theta)**2 + 1./pellet_ellipse**2 * sin(pellet_theta)**2)*(R-pellet_R)**2  &
                + (sin(pellet_theta)**2 + 1./pellet_ellipse**2 * cos(pellet_theta)**2)*(Z-pellet_Z)**2  &
                + 2.*(R-pellet_R)*(Z-pellet_Z)*sin(pellet_theta)*cos(pellet_theta) * (1./pellet_ellipse**2 - 1.) )

  atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
  atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))
  atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))

!  pellet_volume = PI * pellet_radius**2 * pellet_R * pellet_phi ! simulated pellet volume

  phys_pellet_volume = pellet_particles /pellet_density         ! physical pellet volume 

! the number of particles ablated from the physical pellet in units 10^20 m^-3 per unit of JOREK time

  pellet_atomic = 2.d0

!----------------- model from Kuteev (see Polevoi, PPCF2008)
 ! ablation_rate = 1.62d5 * central_density**(-0.77) * T0**(1.72) * r0**(0.45) * phys_pellet_volume**(0.48) * pellet_atomic**0.217
!----------------- NGS model Parks (see Gal NF2008)
  ablation_rate = 2.01d4 * central_density**(-0.81) * max(T0,0.d0)**(1.64) * max(r0,0.d0)**(0.33) * phys_pellet_volume**(0.44) * pellet_atomic**0.5

! particle source in JOREK normalisation

  particle_source = ablation_rate / central_density  * atn * atn_phi / pellet_volume

  volume_source   = atn * atn_phi

  if(volume_source .lt. 0.d0) then
!    print*, 'volume_source is negative. volume_source=', volume_source
  endif

  if(ablation_rate .lt. 0.d0) then
!    print*, 'ablation rate is negative. step.'
  endif
  

end if

return
end subroutine pellet_source2

subroutine update_pellet(my_id,node_List,element_list)
!******************************************************************************
! routine updates the pellet position and the size of the simulated           *
! and physical pellet sizes (from the integral of the pellet particle source) *
!******************************************************************************

use constants
use data_structure
use phys_module
use mpi_mod
implicit none


type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: psi_axis, psi_bnd
integer :: my_id, ierr
real*8  :: V_normalisation, density, density_in, density_out, pressure,pressure_in,pressure_out

real*8 :: R_out, Z_out, s_out, t_out, P0_s,P0_t,P0_st,P0_ss,P0_tt
integer :: i_elm, ifail

if (pellet_amplitude .gt. 0) return

call Integrals_3D(my_id, node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO) ! assumes Deuterium!

pellet_R = pellet_R + pellet_velocity_R * tstep / V_normalisation
pellet_Z = pellet_Z + pellet_velocity_Z * tstep / V_normalisation


phys_ablation = total_pellet_particles * central_density / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO)

total_pellet_particles = total_pellet_particles * central_density * tstep 
total_plasma_particles = total_plasma_particles * central_density          ! undo normalisation


call find_RZ(node_list,element_list,pellet_R,pellet_Z,R_out,Z_out,i_elm,s_out,t_out,ifail)
call interp(node_list,element_list,i_elm,1,1,s_out,t_out,pellet_psi,P0_s,P0_t,P0_st,P0_ss,P0_tt)

if (my_id .eq. 0) then

    pellet_particles = max(pellet_particles - total_pellet_particles, 0.d0)

    write(*,'(A,4e14.6)') ' pellet (R,Z) =', pellet_R, pellet_Z,pellet_velocity_R/V_normalisation,pellet_velocity_Z/V_normalisation
    write(*,'(A,6e14.6,A)') ' total particles added in this step : ', pellet_R, pellet_Z,pellet_particles,total_pellet_particles, total_plasma_particles,phys_ablation,' [10^20]'
    write(*,'(A,4e14.6)') ' remaining particles in pellet      : ', pellet_particles
    write(*,'(A,4e14.6)') ' pellet volume (sim,phys)           : ', total_pellet_volume,pellet_particles/pellet_density

else 
  pellet_particles = 0.0
end if
 
call MPI_Bcast(pellet_particles,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

return 
 
end subroutine update_pellet

subroutine update_spi(my_id,node_List,element_list)
!******************************************************************************
! routine updates the shattered pellet position and the size of the simulated *
! and physical pellet sizes (from the integral of the pellet particle source) *
!******************************************************************************

use constants
use data_structure
use phys_module
use mpi_mod
use mgi_module
implicit none


type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: psi_axis, psi_bnd
integer :: my_id, ierr
real*8  :: V_normalisation, density, density_in, density_out, pressure,pressure_in,pressure_out

real*8 :: R_out, Z_out, s_out, t_out, P0_s,P0_t,P0_st,P0_ss,P0_tt
integer :: i_elm, ifail

V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO) ! assumes Deuterium!

!spi_R = spi_R + spi_Vel_R * tstep / V_normalisation
!spi_Z = spi_Z + spi_Vel_Z * tstep / V_normalisation

     spi_R = mgi_R + (t_now-t_mgi)*spi_Vel_R/V_normalisation
     spi_Z = mgi_Z + (t_now-t_mgi)*spi_Vel_Z/V_normalisation

    write(*,'(A,4e14.6)') ' pellet (R,Z) =', spi_R, spi_Z,spi_Vel_R/V_normalisation,spi_Vel_Z/V_normalisation


return 
 
end subroutine update_spi

end module pellet_module

!> -------------------------------------------------------------------------------------------!
!> The jorek2_particles_radiation post-processor generates synthetic images 
!> and/or radiation signals from jorek2_particles simulations.
!> The assumption used in the jorek2_particles_radiation are:
!> 	-> geometric optics
!>	-> point light sources
!>	-> no occulding objects and no scattering events (direct emission only)
!> The camera models currently implemented are:
!>	-> ideal pinhole camera: total flux and average radiance importance
!> The lightsource models currently implemented are:
!>	-> runaway electron synchrotron emission [J. Schwinger,Phys.Rev.,vol.75,p.1912,1949]
!> -------------------------------------------------------------------------------------------!
program jorek2_particles_radiation
  use mod_particle_sim 
  implicit none
  contains

  ! Variable declarations -------------!
  class(particle_sim) ::   kinetic_sim             !< the particle simulation
  integer ::               num_sims,num_kin_groups !< total number of simulations and groups
  ! -----------------------------------!

  ! Hard-coded values -----------------!
  ! -----------------------------------!

  ! Initialise ------------------------!
  ! -----------------------------------!

  ! Finalise --------------------------!
  call kinetic_sim%finalize()
  ! -----------------------------------!

end program jorek2_particles_radiation

!> Compute and store the ligth soutce properties ---------------------------------------------!
!> inputs:
!>   kinetic_sim:     (particle_sim) the kinetic simulation object 
!>   num_sims:        (int) number of simulations to be loaded
!>   num_kin_groups:  (int) number of kinetic groups of the simulation
!>   filenames:       (array2D(char)) array 2D containing the kinetic sims filename to be loaded
!>   skip_jorek2help: (logical,optional) if present and .false. jorek2help is displayed
!> outputs: 
!>   kinetic_sims: (array(particle_sim)) array containing the kinetic simulations
subroutine initialise_lights(kinetic_sim,num_sims,num_kin_groups,filenames,skip_jorek2help)
  use mod_particle_sim
  use mod_particle_io, only: read_simulation_hdf5,get_simulation_hdf5_time
  implicit none
  contains

  ! Variable declarations -------------!
  class(particle_sim),intent(inout) ::     kinetic_sims
  integer,intent(in) ::                    num_sims
  integer,intent(in) ::                    num_kin_groups
  character,allocatable(:,:),intent(in) :: filenames
  logical,optional,intent(in) ::           skip_jorek2help
  integer ::                               i,j,k
  ! -----------------------------------!

  ! Initialisations -------------------!
  call kinetic_sim%initialize(num_kin_groups, skip_jorek2help)
  ! -----------------------------------!

  ! loop on the simulations
  loop_over_sims:do i=1,num_sims
    call read_simulation_hdf5(kinetic_sim,trim(filenames(:,i))) 
  end do loop_over_sims
  

end subroutine read_kinetic_sims
!> -------------------------------------------------------------------------------------------!

!> Module containing a datatype for simulation parameters and routines for
!> reading this from an input namelist
module mod_particle_sim
  type type_particle_sim
    integer :: i = 1           !< Current step number of the simulation
    real*8  :: t = 0.d0        !< Current time
    real*8  :: dt              !< Timestep

    integer :: iend            !< End at step iend
    real*8  :: tend            !< End at time tend
  end type type_particle_sim

  type type_particle_init_params
    character(len=80) :: particle_restart_file = '' !< Particle restart file to read from

    integer :: species = 0     !< Atomic number Z of the particles (-1) for electrons
    real*4  :: atomic_mass     !< Atomic mass in a.m.u.

    character(len=80) :: location_accept_function = 'location_accept_any' !< Which function to use for particle position rejection sampling
    real*4  :: location_accept_parameters(1:9) = 0 !< Extra arguments for this function
    integer :: particle_seed = 0 !< Seed for PCG random sequence used for particle init
    character(len=6) :: particle_initializer = 'pcg32' !< Method to use for seeding particles (options: pcg32, sobol)
  end type type_particle_init_params


  type type_ion_rec_params
    logical :: particle_ion_rec = .true. !< Calculate ionisation/recombination events
    character(len=6)  :: adas_suffix = '' !< Suffix for adas files to read in (ex: scd50_w.dat => 50_w)
  end type type_ion_rec_params

  type type_jorek_restart_file_params
    integer :: t_particles_begin = -1 !< Number of first JOREK restart file (if -1 use only jorek_restart)
    integer :: t_particles_end  !< Number of last JOREK restart file
  end type type_jorek_restart_file_params
end module mod_particle_sim

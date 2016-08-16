!> A group containing some particles in allocatable data storage.
!> This could be improved in the future by using a different data structure for storage
!> A group has parameters for initialisation and a variable number of pusher hooks
module mod_particle_group
  use mod_particle_base
  implicit none

  type particle_init_params
    integer*1 :: species = 0     !< Atomic number Z of the particles (-1) for electrons
    real*4    :: atomic_mass     !< Atomic mass in a.m.u.

    character(len=80) :: location_accept_function = 'location_accept_any' !< Which function to use for particle position rejection sampling
    real*4  :: location_accept_parameters(1:9) = 0 !< Extra arguments for this function
    integer :: particle_seed = 0 !< Seed for PCG random sequence used for particle init
    character(len=6) :: particle_initializer = 'pcg32' !< Method to use for seeding particles (options: pcg32, sobol)
  end type particle_init_params


  type ion_rec_params
    character(len=6)  :: adas_suffix = '' !< Suffix for adas files to read in (ex: scd50_w.dat => 50_w)
  end type ion_rec_params





  type :: particle_group
    class(particle_base), dimension(:), allocatable :: particles
    type(particle_init_params)                      :: init_params
  end type particle_group


end module mod_particle_group

!> A group containing some particles in allocatable data storage.
!> This could be improved in the future by using a different data structure for storage
!> A group has parameters for initialisation and a variable number of pusher hooks
module mod_particle_group
  use mod_particle_base
  implicit none

  type ion_rec_params
    character(len=6)  :: adas_suffix = '' !< Suffix for adas files to read in (ex: scd50_w.dat => 50_w)
  end type ion_rec_params

  type :: particle_group
    class(particle_base), dimension(:), allocatable :: particles
  end type particle_group

end module mod_particle_group

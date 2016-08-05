program particles

use mod_boris ! or mod_gc_NAME or something similar
implicit none

type(particle_boris), allocatable, dimension(:) :: particles
type(particle_sim) :: sim
class(fields), allocatable :: fields
class(particle_hdf5_io), allocatable :: io

integer :: i


call setup

! Prepare HDF5 IO
allocate(particle_hdf5_io_boris::io)
allocate(particle_boris::io%particle_type)
call io%set_data_type

call read_parameters(sim)
call allocate(particle_boris::particles(sim%n_particles))
call read_adas
call initialize_rng(rng)
call initialize_particles(sim)

!$omp parallel
do sim%i=1,iend
  !$omp do
  do i=1,size(particles,1)
    call push_boris(fields, particles(i), sim)
  enddo
  !$end omp do

  !$omp single
  call diagnostics
  !$end omp single
  sim%t = sim%t + sim%dt
enddo
!$end omp parallel

call teardown

end program particles

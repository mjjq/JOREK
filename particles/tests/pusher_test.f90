!> Program to test different pushers on computation time and accuracy
!> Be sure to set OMP_NUM_THREADS to the number of cores your system has
!> to prevent hyperthreading from messing up results, or disable openmp completely
!> for the best results.
!> Requires gnuplot to generate figures
program pusher_test
  use mod_penning_case
  use mod_gradB_case
  use mod_constants, only: CARTESIAN, CYLINDRICAL
  implicit none
  real*8, dimension(4), parameter :: penning_timesteps = [1d-2, 1d-3, 1d-4, 1d-5]
  real*8, dimension(4), parameter :: gradB_timesteps = [1d-2, 1d-3, 1d-4, 1d-5]

  call boris_cases

contains




subroutine boris_cases
  use mod_boris
  type(case_penning)   :: penning_cartesian
  type(case_penning)   :: penning_cylindrical
  type(case_gradB)     :: gradB_cartesian
  type(case_gradB)     :: gradB_cylindrical
  type(particle_boris) :: particle
  type(pusher_boris)   :: pusher
  integer :: i
  real*8 :: err, runtime
  
  ! Setup the cases (sets the fields)
  penning_cartesian   = case_penning(geometry=CARTESIAN)
  penning_cylindrical = case_penning(geometry=CYLINDRICAL)
  gradB_cartesian     = case_gradB(geometry=CARTESIAN)
  gradB_cylindrical   = case_gradB(geometry=CYLINDRICAL)

  !$omp parallel default(shared) &
  !$omp private(i, err, runtime, pusher, particle)
  !$omp do
  do i=size(penning_timesteps),1,-1
    pusher%fixed_timestep = penning_timesteps(i)
    call penning_cartesian%initialize_particle(particle, pusher)
    call penning_cartesian%run(pusher, particle, err, runtime)
    write(*,*) "Penning Boris_Cartesian  ", pusher%fixed_timestep, err, runtime
  end do
  !$omp end do
  
  !$omp do
  do i=size(penning_timesteps),1,-1
    pusher%fixed_timestep = penning_timesteps(i)
    call penning_cylindrical%initialize_particle(particle, pusher)
    call penning_cylindrical%run(pusher, particle, err, runtime)
    write(*,*) "Penning Boris_Cylindrical", pusher%fixed_timestep, err, runtime
  end do
  !$omp end do

  !$omp do
  do i=size(gradB_timesteps),1,-1
    pusher%fixed_timestep = gradB_timesteps(i)
    call gradB_cartesian%initialize_particle(particle, pusher)
    call gradB_cartesian%run(pusher, particle, err, runtime)
    write(*,*) "gradB Boris_Cartesian", pusher%fixed_timestep, err, runtime
  end do
  !$omp end do

  !$omp do
  do i=size(gradB_timesteps),1,-1
    pusher%fixed_timestep = gradB_timesteps(i)
    call gradB_cylindrical%initialize_particle(particle, pusher)
    call gradB_cylindrical%run(pusher, particle, err, runtime)
    write(*,*) "gradB Boris_Cylindrical", pusher%fixed_timestep, err, runtime
  end do
  !$omp end do
  !$omp end parallel

  ! Run a single case with output
  !pusher%fixed_timestep = 1d-3
  !call gradB_cartesian%initialize_particle(particle, pusher)
  !call gradB_cartesian%run(pusher, particle, err, runtime, output_file='gradB_cartesian_1d-3.txt')
end subroutine boris_cases
end program pusher_test

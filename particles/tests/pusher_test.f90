!> Program to test different pushers on computation time and accuracy
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
  integer :: i, u, stat
  real*8 :: err, runtime
  
  ! Setup the cases (sets the fields)
  penning_cartesian   = case_penning(geometry=CARTESIAN)
  penning_cylindrical = case_penning(geometry=CYLINDRICAL)
  gradB_cartesian     = case_gradB(geometry=CARTESIAN)
  gradB_cylindrical   = case_gradB(geometry=CYLINDRICAL)

  open(newunit=u, file="boris.txt")
  do i=1,size(penning_timesteps)
    pusher%fixed_timestep = penning_timesteps(i)
    call penning_cartesian%initialize_particle(particle, pusher)
    call penning_cartesian%run(pusher, particle, err, runtime)
    write(u,*) "Penning Boris_Cartesian", pusher%fixed_timestep, err, runtime
  end do
  do i=1,size(penning_timesteps)
    pusher%fixed_timestep = penning_timesteps(i)
    call penning_cylindrical%initialize_particle(particle, pusher)
    call penning_cylindrical%run(pusher, particle, err, runtime)
    write(u,*) "Penning Boris_Cylindrical", pusher%fixed_timestep, err, runtime
  end do
  do i=1,size(gradB_timesteps)
    pusher%fixed_timestep = gradB_timesteps(i)
    call gradB_cartesian%initialize_particle(particle, pusher)
    call gradB_cartesian%run(pusher, particle, err, runtime)
    write(u,*) "gradB Boris_Cartesian", pusher%fixed_timestep, err, runtime
  end do
  do i=1,size(gradB_timesteps)
    pusher%fixed_timestep = gradB_timesteps(i)
    call gradB_cylindrical%initialize_particle(particle, pusher)
    call gradB_cylindrical%run(pusher, particle, err, runtime)
    write(u,*) "gradB Boris_Cylindrical", pusher%fixed_timestep, err, runtime
  end do
  close(u)
  ! Plot the above results
  call system('gnuplot -e "set logscale xy; set xlabel \"timestep size\"; &
      set terminal png; &
      set key bottom right; &
      set format x \"10^{%L}\"; set format y \"10^{%L}\"; &
      set output \"media/tests/all_pushers/penning.png\"; &
      plot \"< grep Penning boris.txt | grep Cartesian\" u 3:4 w l t \"Boris Cartesian\", &
           \"< grep Penning boris.txt | grep Cylindrical\" u 3:4 w l t \"Boris Cylindrical\"; &
      set output \"media/tests/all_pushers/gradB.png\"; &
      plot \"< grep gradB boris.txt | grep Cartesian\" u 3:4 w l t \"Boris Cartesian\", &
           \"< grep gradB boris.txt | grep Cylindrical\" u 3:4 w l t \"Boris Cylindrical\" &
      "')

  ! Run a single case with output (for demo) and plot the result
  pusher%fixed_timestep = 1d-2
  call gradB_cartesian%initialize_particle(particle, pusher)
  call gradB_cartesian%run(pusher, particle, err, runtime, output_file='gradB_cartesian_1d-2.txt')
  call system('gnuplot -e "set xlabel \"x\"; set ylabel \"y\"; &
      set terminal png; &
      set output \"media/tests/gradB/gradB_xy_boris.png\"; &
      plot \"gradB_cartesian_1d-2.txt\" u 1:2 w l"')

  ! Run a single case with output (for demo) and plot the result
  pusher%fixed_timestep = 1d-2
  call penning_cartesian%initialize_particle(particle, pusher)
  call penning_cartesian%run(pusher, particle, err, runtime, output_file='penning_cartesian_1d-2.txt')
  call system('gnuplot -e "set xlabel \"x\"; set ylabel \"y\"; &
      set terminal png; &
      set output \"media/tests/penning/penning_xy_boris.png\"; &
      plot \"penning_cartesian_1d-2.txt\" u 1:2 w l"')

  ! delete the created files again
  open(newunit=u, iostat=stat, file='boris.txt', status='old')
  if (stat .eq. 0) close(u, status='delete')
  open(newunit=u, iostat=stat, file='gradB_cartesian_1d-2.txt', status='old')
  if (stat .eq. 0) close(u, status='delete')
  open(newunit=u, iostat=stat, file='penning_cartesian_1d-2.txt', status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine boris_cases
end program pusher_test

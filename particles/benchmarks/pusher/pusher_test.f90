!> Program to test different pushers on computation time and accuracy.
!> Requires gnuplot to generate figures.
!> Run with `make docs`
!>
!># Results
!>## Penning test (see [[mod_penning_case]])
!> ![pusher-test-penning](|media|/tests/all_pushers/penning.png)
!> ![pusher-test-penning-time](|media|/tests/all_pushers/penning_time.png)
!>## gradB test (see [[mod_gradB_case]])
!> ![pusher-test-gradB](|media|/tests/all_pushers/gradB.png)
program pusher_test
  use mod_particle_types
  use mod_penning_case
  use mod_gradB_case
  implicit none
  real*8, dimension(4), parameter :: penning_timesteps = [1d-2, 1d-3, 1d-4, 1d-5]
  real*8, dimension(4), parameter :: gradB_timesteps = [1d-2, 1d-3, 1d-4, 1d-5]

  call test_boris
contains


!> Run some tests with the boris pusher.
subroutine test_boris
  use mod_boris
  type(particle_kinetic_leapfrog) :: particle
  integer :: i, u, stat, k, n_steps
  real*8 :: err, t0, t1
  type(case_penning_cylindrical) :: case1
  type(case_penning_cartesian)   :: case2
  type(case_gradB_cylindrical)   :: case3
  type(case_gradB_cartesian)     :: case4
  character(len=*), parameter :: format = "(A,g9.2,g16.8,g16.8)"

  open(newunit=u, file="boris.txt")
  do i=1,size(penning_timesteps)
    associate (c => case1, dt => penning_timesteps(i))
      n_steps = nint(c%time_end/dt)
      call c%initialize(particle)
      call boris_initial_half_step_backwards(particle, c%E(particle%x, 0.d0), &
        c%B(particle%x, 0.d0), dt)
      call cpu_time(t0)
      do k=1,n_steps
        call boris_push_cylindrical(particle, c%E(particle%x, 0.d0), c%B(particle%x, 0.d0), dt)
      end do
      call cpu_time(t1)
      err = c%calc_error(particle)
      write(u,format) "Penning Boris_Cylindrical", dt, err, t1-t0
    end associate
  end do
  do i=1,size(penning_timesteps)
    associate (c => case2, dt => penning_timesteps(i))
      n_steps = nint(c%time_end/dt)
      call c%initialize(particle)
      call boris_initial_half_step_backwards(particle, c%E(particle%x, 0.d0), &
        c%B(particle%x, 0.d0), dt)
      call cpu_time(t0)
      do k=1,n_steps
        call boris_push_cartesian(particle, c%E(particle%x, 0.d0), c%B(particle%x, 0.d0), dt)
      end do
      call cpu_time(t1)
      err = c%calc_error(particle)
      write(u,format) "Penning Boris_Cartesian", dt, err, t1-t0
    end associate
  end do
  do i=1,size(gradB_timesteps)
    associate (c => case3, dt => gradB_timesteps(i))
      n_steps = nint(c%time_end/dt)
      call c%initialize(particle)
      call boris_initial_half_step_backwards(particle, c%E(particle%x, 0.d0), &
        c%B(particle%x, 0.d0), dt)
      call cpu_time(t0)
      do k=1,n_steps
        call boris_push_cylindrical(particle, c%E(particle%x, 0.d0), c%B(particle%x, 0.d0), dt)
      end do
      call cpu_time(t1)
      err = c%calc_error(particle)
      write(u,format) "gradB Boris_Cylindrical", dt, err, t1-t0
    end associate
  end do
  do i=1,size(gradB_timesteps)
    associate (c => case4, dt => gradB_timesteps(i))
      n_steps = nint(c%time_end/dt)
      call c%initialize(particle)
      call boris_initial_half_step_backwards(particle, c%E(particle%x, 0.d0), &
        c%B(particle%x, 0.d0), dt)
      call cpu_time(t0)
      do k=1,n_steps
        call boris_push_cartesian(particle, c%E(particle%x, 0.d0), c%B(particle%x, 0.d0), dt)
      end do
      call cpu_time(t1)
      err = c%calc_error(particle)
      write(u,format) "gradB Boris_Cartesian", dt, err, t1-t0
    end associate
  end do
  close(u)

  ! Run a single case with output (for demo) and write the results
  open(newunit=u, file='gradB_cartesian_1d-2.txt')
  associate (c => case4, dt => 1d-2)
    n_steps = nint(c%time_end/dt)
    call c%initialize(particle)
    call boris_initial_half_step_backwards(particle, c%E(particle%x, 0.d0), &
      c%B(particle%x, 0.d0), dt)
    call cpu_time(t0)
    do k=1,n_steps
      call boris_push_cartesian(particle, c%E(particle%x, 0.d0), c%B(particle%x, 0.d0), dt)
      if (mod(k,10) .eq. 0) write(u,*) particle%x
    end do
    call cpu_time(t1)
    err = c%calc_error(particle)
  end associate
  close(u)
  ! Run a single case with output (for demo) and write the results
  open(newunit=u, file='penning_cartesian_1d-2.txt')
  associate (c => case2, dt => 1d-2)
    n_steps = nint(c%time_end/dt)
    call c%initialize(particle)
    call boris_initial_half_step_backwards(particle, c%E(particle%x, 0.d0), &
      c%B(particle%x, 0.d0), dt)
    call cpu_time(t0)
    do k=1,n_steps
      call boris_push_cartesian(particle, c%E(particle%x, 0.d0), c%B(particle%x, 0.d0), dt)
      if (mod(k,10) .eq. 0) write(u,*) particle%x
    end do
    call cpu_time(t1)
    err = c%calc_error(particle)
  end associate
  close(u)

  ! Create a temporary file for gnuplot
  open(newunit=u, file='pusher_test.gp', status='new')
  write(u,"(A)") '&
      set terminal png; &
      set key bottom right; &
      set logscale xy; set xlabel "timestep size"; &
      set format x "10^{%L}"; set format y "10^{%L}"; &
      set output "media/tests/all_pushers/penning.png"; &
      set title "Penning trap test"; &
      plot "< grep Penning boris.txt | grep Cartesian" u 3:4 w l t "Boris Cartesian", &
           "< grep Penning boris.txt | grep Cylindrical" u 3:4 w l t "Boris Cylindrical"; &
      set ylabel "calculation time [s]"; &
      set output "media/tests/all_pushers/penning_time.png"; &
      plot "< grep Penning boris.txt | grep Cartesian" u 5:4 w l t "Boris Cartesian", &
           "< grep Penning boris.txt | grep Cylindrical" u 5:4 w l t "Boris Cylindrical"; &
      set title "gradB test"; &
      set xlabel "timestep size"; &
      set output "media/tests/all_pushers/gradB.png"; &
      plot "< grep gradB boris.txt | grep Cartesian" u 3:4 w l t "Boris Cartesian", &
           "< grep gradB boris.txt | grep Cylindrical" u 3:4 w l t "Boris Cylindrical"; &

      set xlabel "x"; set ylabel "y"; &
      unset logscale; unset format x; unset format y; &
      set output "media/tests/gradB/gradB_xy_boris.png"; &
      plot "gradB_cartesian_1d-2.txt" u 1:2 w l; &
      set output "media/tests/penning/penning_xy_boris.png"; &
      plot "penning_cartesian_1d-2.txt" u 1:2 w l'
  close(u)

  call system('gnuplot pusher_test.gp')

  ! delete the created files again
  open(newunit=u, iostat=stat, file='pusher_test.gp', status='old')
  if (stat .eq. 0) close(u, status='delete')
  open(newunit=u, iostat=stat, file='boris.txt', status='old')
  if (stat .eq. 0) close(u, status='delete')
  open(newunit=u, iostat=stat, file='gradB_cartesian_1d-2.txt', status='old')
  if (stat .eq. 0) close(u, status='delete')
  open(newunit=u, iostat=stat, file='penning_cartesian_1d-2.txt', status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine test_boris
end program pusher_test

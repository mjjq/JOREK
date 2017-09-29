!> This program calculates the charge state equilibrium with a MC method
!> and with the coronal equilibrium, and compares the two.
!> Requires gnuplot to generate figures.
!> Run with `make docs`
!>
!># Results
!>## Charge state distribution with temperature
!> ![charge-state-temperature](|media|/tests/openadas/charge_state_temperature.png)
!>## Equilibration in time (MC method) and comparison to coronal equilibrium value
!> ![charge-state-time](|media|/tests/openadas/charge_state_time.png)
!>
!> We compare here many timesteps of the coronal equilibrium calculation
!> but the results are very similar to those obtained by using just one step.
program compare_mc_coronal
use mod_openadas
use mod_coronal
use mod_ionisation_recombination
use mpi
use mod_pcg32_rng
use mod_random_seed
implicit none

character(len=6), parameter :: suffix = "50_w"
character(len=80), parameter :: dir = "particles/examples/"
type(ADF11_all) :: ad
integer :: ierr, provided
call MPI_INIT_Thread(MPI_THREAD_SINGLE, provided, ierr)

ad = read_adf11(suffix, directory=dir)

call plot_coronal_equilibrium(ad)
call plot_mc_coronal_equilibria(ad)

call MPI_Finalize(ierr)
contains

subroutine plot_mc_coronal_equilibria(ad)
  type(ADF11_all), intent(in)  :: ad
  character(len=80), parameter :: datafile_mc = "W_mc_equilibrium.txt"
  character(len=80), parameter :: datafile_cor = "W_coronal_equilibria.txt"
  character(len=80), parameter :: plotfile = "plot_mc_coronal_equilibria.gp"

  real*8, parameter  :: tstep = 1.d-6 !< [s]
  integer, parameter :: nstep = 4000
  integer, parameter :: n_particles = 1000
  real*8, parameter  :: density = 20.d0 !< log10 density in m^-3
  real*8, parameter  :: temperature = 7.677 !< log10 temperature in K: 7.677 = 4.1 keV

  integer, dimension(:), allocatable :: z !< Charge states of each particle
  integer, dimension(:), allocatable :: p_mc
  real*8, dimension(:), allocatable  :: p_cor !< population at each charge state at this time
  integer :: i, it, u_mc, u_cor, u

  type(pcg32_rng) :: rng
  real*8 :: ran2(2)
  call rng%initialize(2, random_seed(), 1, 1)

  allocate(p_mc(0:ad%n_Z))
  allocate(p_cor(0:ad%n_Z))
  allocate(z(n_particles))

  open(newunit=u_mc,  file=datafile_mc,  status="replace")
  open(newunit=u_cor, file=datafile_cor, status="replace")

  ! Set z over all particles from 1 to adf11%n_Z
  do i=1,n_particles
    z(i) = nint(ad%n_Z * real(i,8)/real(n_particles,8))
  enddo

  p_cor = 1.d0
  ! Write output in format:
  ! ROW1x COL1y VAL
  ! (blank)
  ! ROW2x COL1y VAL
  do i=1,nstep
    p_mc = 0
    do it=1,n_particles
      call rng%next(ran2)
      z(it) = new_charge(z(it), ad, density, temperature, tstep, ran2)
      p_mc(z(it)) = p_mc(z(it)) + 1
    enddo
    do it=lbound(p_mc,1),ubound(p_mc,1)
      write(u_mc,"(3g16.8)") i*tstep, it, real(p_mc(it),8)/real(sum(p_mc),8)
    enddo

    call coronal_timestep(ad, p_cor, tstep, density, temperature)
    do it=lbound(p_cor,1),ubound(p_cor,1)
      write(u_cor,"(3g16.8)") i*tstep, it, p_cor(it)/sum(p_cor)
    enddo
  enddo
  close(u_mc)
  close(u_cor)

  ! Plot the results
  open(newunit=u, file=plotfile, status='replace')
  write(u,"(A)") '&
      set xr [0:74]; &
      set yr [0:4e-3]; &
      set cbrange [0:0.4]; &
      set terminal png; &
      set key top right; &
      set xlabel "Z"; &
      set output "media/tests/openadas/charge_state_time.png"; &
      set multiplot layout 2,1'
  write(u,"(A,A,A)") 'plot "', trim(datafile_mc), '" using 2:1:3 with image t "MC"'
  write(u,"(A,A,A)") 'plot "', trim(datafile_cor), '" using 2:1:3 with image t "COR"'
  close(u)

  call system('gnuplot '//plotfile)

  call rm(datafile_mc)
  call rm(datafile_cor)
  call rm(plotfile)
end subroutine

subroutine plot_coronal_equilibrium(ad)
  type(ADF11_all), intent(in)    :: ad
  real*8, dimension(0:ad%n_Z)    :: p
  real*8 :: T
  integer :: t_exp, u, i
  character(len=80), parameter :: datafile = "W_coronal_equilibrium.txt"
  character(len=80), parameter :: plotfile = "plot_coronal_equilibrium.gp"
  character(len=3) :: Z

  open(newunit=u, file=datafile, status="replace")
  do t_exp = 1000,8000
    T = real(t_exp,8)/1000.d0 ! in log10 K

    p = 1.d0 ! initialize to 1
    call coronal_timestep(ad, p, 1.d0, 20.d0, T) ! use a fixed large timestep of 1 to solve
    write(u,"(100g16.8)") 10.d0**T, p/sum(p)
  end do
  close(u)

  ! Plot the results
  open(newunit=u, file=plotfile, status='replace')
  write(u,"(A)") '&
      set terminal png; &
      set key top right; &
      set logscale xy; set xlabel "T [K]"; &
      set format x "10^{%L}"; set format y "10^{%L}"; &
      set yrange [0.0000001:1]; &
      set output "media/tests/openadas/charge_state_temperature.png"; &
      plot \'
  do i=1,ad%n_Z
    write(Z, "(i3)") i
    write(u, "(5A)") '"', trim(datafile), '" using 1:', adjustl(Z), ' with lines notitle, \'
  end do
  close(u)

  call system('gnuplot '//plotfile)

  ! Delete the temporary file
  call rm(datafile)
  call rm(plotfile)
end subroutine plot_coronal_equilibrium

subroutine rm(filename)
  character(len=*), intent(in) :: filename
  integer :: u, stat
  open(newunit=u, iostat=stat, file=filename, status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine rm
end program compare_mc_coronal

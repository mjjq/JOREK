!> This module contains some testcases for projecting particles.
!> First we test the projection and scaling with n.
!> Then we look into the effect of the smoothing parameter s.
module particle_projection_spec
use mod_project_particles
use data_structure
use mod_particle_types
use mod_pcg32_rng
use mod_sobseq_rng
use projection_spec, only: elements_mean_rms, f_1, default_square_grid, default_polar_grid, &
    default_flux_grid_32, default_flux_grid_31, calc_rhs_f
use constants, only: TWOPI
use fruit
implicit none

logical, parameter :: write_proj_output = .false. !< Set to true to write restart files with the projected density
logical, parameter :: EXTRATEST = .false. !< Set to .true. to do flux-aligned grid projection tests

contains

subroutine setup_particle_projection_spec
  call initialise_basis
end subroutine setup_particle_projection_spec

!> Project 10^3-10^5 particles generated with pcg onto square grid
subroutine test_square_10_10_pcg
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_n_square_10_10(node_list, element_list, [1000,10000,100000], pcg32_rng(), mean_tol=3d-8*[1,1,1], rms_tol=23d0/sqrt([1d3,1d4,1d5]))
end subroutine test_square_10_10_pcg

!> Project 10^3-10^5 particles generated with sobseq onto square grid
subroutine test_square_10_10_sob
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_n_square_10_10(node_list, element_list, [1000,10000,100000], sobseq_rng(), mean_tol=3d-8*[1,1,1], rms_tol=405d0/[1d3,1d4,1d5])
end subroutine test_square_10_10_sob


!> Project 10^3-10^5 particles generated with sobseq onto polar grid
subroutine test_polar_30_22_sob
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_polar_grid(node_list, element_list, 22)
  call project_n_polar_30_npol(node_list, element_list, [1000,10000,100000], sobseq_rng(), mean_tol=3d-5*[1,1,1], rms_tol=5d4/[1d3,1d4,1d5])
end subroutine test_polar_30_22_sob
!> Project 10^3-10^5 particles generated with sobseq onto polar grid
subroutine test_polar_30_21_sob
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_polar_grid(node_list, element_list, 21)
  call project_n_polar_30_npol(node_list, element_list, [1000,10000,100000], sobseq_rng(), mean_tol=3d-5*[1,1,1], rms_tol=5d4/[1d3,1d4,1d5])
end subroutine test_polar_30_21_sob

!> Project 10^3-10^5 particles generated with sobseq onto flux grid (odd)
subroutine test_flux_40_31_pcg
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  real*8 :: w
  if (.not. EXTRATEST) return
  call default_flux_grid_31(node_list, element_list)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [1000,10000,100000], pcg32_rng(), 'flux_40_31', volume=w, mean_tol=2d-5*[1,1,1], rms_tol=45d0/sqrt([1d3,1d4,1d5]))
end subroutine test_flux_40_31_pcg
!> Project 10^3-10^5 particles generated with sobseq onto flux grid (even)
subroutine test_flux_40_32_pcg
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  real*8 :: w
  if (.not. EXTRATEST) return
  call default_flux_grid_32(node_list, element_list)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [1000,10000,100000], pcg32_rng(), 'flux_40_32', volume=w, mean_tol=2d-5*[1,1,1], rms_tol=45d0/sqrt([1d3,1d4,1d5]))
end subroutine test_flux_40_32_pcg


!> Create a RHS by integrating f_1 and with monte carlo methods and check that they are close
!> This guards against errors in node indices etc
subroutine rhs_convergence_square_10_10(node_list, element_list, n, p, rng)
  use mod_initialise_particles
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout)  :: element_list
  type(DMUMPS_STRUC), intent(inout) :: p
  integer, intent(in) :: n
  class(type_rng), intent(in) :: rng

  integer :: n_AA, ielm_out, ifail, i
  real*8 :: R_out, Z_out, s_out, t_out, tol
  real*8, allocatable, dimension(:) :: rhs_f
  type(particle_fieldline), dimension(:), allocatable :: particles
  character*8 :: n_s
  write(n_s, '(i8)') n

  n_AA = maxval(node_list%node(1:node_list%n_nodes)%index(4))
  allocate(rhs_f(n_AA))
  allocate(particles(n))

  ! to prevent omp trouble
  call find_RZ(node_list,element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)
  call initialise_particles(particles, node_list, element_list, rng)
  particles(:)%weight = TWOPI/real(n)
  call project_particles(node_list, element_list, p, particles, 1, skip_proj=.true.)
  call calc_rhs_f(node_list,element_list,f_1,rhs_f)
  call assert_false(isnan(sum(rhs_f)), 'sum integrated rhs is not nan[n='//trim(adjustl(n_s))//']')
  call assert_false(isnan(sum(p%rhs)), 'sum MC rhs is not nan[n='//trim(adjustl(n_s))//']')
  select type (rng)
  type is (pcg32_rng)
    tol = 7d0/sqrt(real(n))
  type is (sobseq_rng)
    tol = 65d0/real(n)
  end select
  call assert_equals(0.d0, sum(abs(rhs_f-p%rhs)), tol, 'sum abs integrated - MC rhs [n='//trim(adjustl(n_s))//']')
  call assert_equals(0.d0, maxval(abs(rhs_f-p%rhs)), tol, 'max abs integrated - MC rhs [n='//trim(adjustl(n_s))//']')
end subroutine rhs_convergence_square_10_10

!> Test convergence of RHS for n particles
subroutine test_rhs_square_10_10_pcg
  type(DMUMPS_STRUC) :: p
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call prepare_mumps_par(node_list, element_list, p, smoothing=0.d0)
  call rhs_convergence_square_10_10(node_list, element_list, 1000, p, pcg32_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 10000, p, pcg32_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 100000, p, pcg32_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 1000000, p, pcg32_rng())
  p%JOB=-2
  call DMUMPS(p)
end subroutine test_rhs_square_10_10_pcg
!> Test convergence of RHS for n particles
subroutine test_rhs_square_10_10_sob
  type(DMUMPS_STRUC) :: p
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call prepare_mumps_par(node_list, element_list, p, smoothing=0.d0)
  call rhs_convergence_square_10_10(node_list, element_list, 1000, p, sobseq_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 10000, p, sobseq_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 100000, p, sobseq_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 1000000, p, sobseq_rng())
  p%JOB=-2
  call DMUMPS(p)
end subroutine test_rhs_square_10_10_sob


!*************************
! Smoothing parameter test
!*************************

!> Test convergence of RHS for 10000 particles with varying smoothing factor
subroutine test_polar_30_22_10000_sob_smoothing
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  real*8 :: w, s, x
  integer :: i
  character(len=8) :: ss
  call default_polar_grid(node_list, element_list, 22)
  w=TWOPI**2*R_geo*amin**2/2.d0
  do i=1,7
    x = real(i-8)
    s = 10d0**x
    write(ss,'(g8.1)') s
    call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s'//ss, volume=w, smoothing=s, &
        ! Calculate error from fit (set rms_tol to 0 to get errors) * 1.2
        rms_tol=[10.d0**(-0.0738*x**2 - 0.972*x - 3.71)*1.2], mean_tol=[2d-5])
  end do
end subroutine test_polar_30_22_10000_sob_smoothing



!> Helper function to project n particles onto a 10x10 square grid
subroutine project_n_square_10_10(node_list, element_list, n, rng, rms_tol, mean_tol)
  use mod_rng
  class(type_rng), intent(in) :: rng
  integer, intent(in), dimension(:) :: n
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, intent(in), dimension(:) :: rms_tol, mean_tol

  call project_n(node_list, element_list, n, rng, 'square_10_10', volume=TWOPI, rms_tol=rms_tol, mean_tol=mean_tol)
end subroutine project_n_square_10_10

!> Helper function to project n particles onto a 30_npol polar grid
subroutine project_n_polar_30_npol(node_list, element_list, n, rng, rms_tol, mean_tol)
  use phys_module
  use mod_rng
  class(type_rng), intent(in) :: rng
  integer, intent(in), dimension(:) :: n
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, intent(in), dimension(:) :: rms_tol, mean_tol

  call project_n(node_list, element_list, n, rng, 'polar_30_32', volume=TWOPI**2*R_geo*amin**2/2.d0, rms_tol=rms_tol, mean_tol=mean_tol)
end subroutine project_n_polar_30_npol

!> Helper function to project n particles onto a grid in node_list, element_list with optional smoothing
subroutine project_n(node_list, element_list, n, rng, name, volume, smoothing, rms_tol, mean_tol)
  use mod_initialise_particles
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  integer, intent(in), dimension(:) :: n
  class(type_rng), intent(in) :: rng
  character(len=*), intent(in) :: name
  real*8, intent(in) :: volume !< total volume
  real*8, intent(in), dimension(:) :: rms_tol, mean_tol
  real*8, optional, intent(in) :: smoothing
  type(particle_kinetic), dimension(:), allocatable :: particles
  type(DMUMPS_STRUC) :: p

  integer :: i, j, ifail, ielm_out
  real*8 :: x(3), R, Z, R_out, Z_out, Phi, tol, s, s_out, t_out
  real*8 :: m, e !< mean, rms error
  character*8 :: n_s, tol_s, ss

  s = 0.d0
  ss = ''
  if (present(smoothing)) then
    s = smoothing
    write(ss,'(g8.1)') s
  end if
  call prepare_mumps_par(node_list, element_list, p, smoothing=s)


  do j=1,size(n)
    write(n_s, '(i8)') n(j)

    allocate(particles(n(j)))
    call find_RZ(node_list,element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)
    call initialise_particles(particles, node_list, element_list, rng)
    particles(:)%weight = volume/real(n(j))

    call project_particles(node_list, element_list, p, particles, 1)

    ! test rms
    call elements_mean_rms(node_list, element_list, f_1, m, e)
    write(tol_s, '(g8.1)') mean_tol(j)
    call assert_equals(1.d0, m, mean_tol(j), 'mean value 1 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']'//trim(ss))
    write(tol_s, '(g8.1)') rms_tol(j)
    call assert_equals(0.d0, e, rms_tol(j), 'rms value 0 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']'//trim(ss))
    if (write_proj_output) then
      call write_particle_distribution_to_h5(node_list, element_list, &
        filename='part_'//trim(adjustl(n_s))//'_'//trim(name)//'.h5', n_fields=1, time=0.d0)
    end if

    deallocate(particles)
  end do
  p%JOB=-2
  call DMUMPS(p)
end subroutine project_n
end module particle_projection_spec

!> This module contains some testcases for projecting particles, ensuring
!> that the projection matrix, RHS and MUMPS work for these cases.
!>
!> we test the projection of a set of particles
!> Note that now the projection seems to converge to the wrong mean, 1% higher than the particle density
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
include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure

logical, parameter :: write_proj_output = .true.
logical, parameter :: EXTRATEST = .false.

contains

subroutine setup_particle_projection_spec
  call initialise_basis
end subroutine setup_particle_projection_spec

subroutine test_square_10_10_pcg
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_n_square_10_10(node_list, element_list, [1000,10000,100000], pcg32_rng())
end subroutine test_square_10_10_pcg

subroutine test_square_10_10_sob
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_n_square_10_10(node_list, element_list, [1000,10000,100000], sobseq_rng())
end subroutine test_square_10_10_sob


subroutine test_polar_30_22_sob
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_polar_grid(node_list, element_list, 22)
  call project_n_polar_30_npol(node_list, element_list, [1000,10000,100000], sobseq_rng())
end subroutine test_polar_30_22_sob
subroutine test_polar_30_21_sob
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_polar_grid(node_list, element_list, 21)
  call project_n_polar_30_npol(node_list, element_list, [1000,10000,100000], sobseq_rng())
end subroutine test_polar_30_21_sob

subroutine test_flux_40_31_pcg
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  real*8 :: w
  if (.not. EXTRATEST) return
  call default_flux_grid_31(node_list, element_list)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [1000,10000,100000], pcg32_rng(), 'flux_40_31', volume=w)
end subroutine test_flux_40_31_pcg
subroutine test_flux_40_32_pcg
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  real*8 :: w
  if (.not. EXTRATEST) return
  call default_flux_grid_32(node_list, element_list)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [1000,10000,100000], pcg32_rng(), 'flux_40_32', volume=w)
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
  real*8 :: R_out, Z_out, s_out, t_out
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
  call assert_equals(0.d0, sum(abs(rhs_f-p%rhs)), 'sum abs integrated - MC rhs [n='//trim(adjustl(n_s))//']')
  call assert_equals(0.d0, maxval(abs(rhs_f-p%rhs)), 'max abs integrated - MC rhs [n='//trim(adjustl(n_s))//']')
end subroutine rhs_convergence_square_10_10

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



subroutine test_polar_30_22_10000_sob_smoothing
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  real*8 :: w
  call default_polar_grid(node_list, element_list, 22)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-7', volume=w, smoothing=1d-7)
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-6', volume=w, smoothing=1d-6)
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-5', volume=w, smoothing=1d-5)
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-4', volume=w, smoothing=1d-4)
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-3', volume=w, smoothing=1d-3)
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-2', volume=w, smoothing=1d-2)
  call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s1d-1', volume=w, smoothing=1d-1)
end subroutine test_polar_30_22_10000_sob_smoothing



subroutine project_n_square_10_10(node_list, element_list, n, rng)
  use mod_rng
  class(type_rng), intent(in) :: rng
  integer, intent(in), dimension(:) :: n
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list

  call project_n(node_list, element_list, n, rng, 'square_10_10', volume=TWOPI)
end subroutine project_n_square_10_10

subroutine project_n_polar_30_npol(node_list, element_list, n, rng)
  use phys_module
  use mod_rng
  class(type_rng), intent(in) :: rng
  integer, intent(in), dimension(:) :: n
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list

  call project_n(node_list, element_list, n, rng, 'polar_30_32', volume=TWOPI**2*R_geo*amin**2/2.d0)
end subroutine project_n_polar_30_npol

subroutine project_n(node_list, element_list, n, rng, name, volume, smoothing)
  use mod_initialise_particles
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  integer, intent(in), dimension(:) :: n
  class(type_rng), intent(in) :: rng
  character(len=*), intent(in) :: name
  real*8, intent(in) :: volume !< total volume
  real*8, optional, intent(in) :: smoothing
  type(particle_kinetic), dimension(:), allocatable :: particles
  type(DMUMPS_STRUC) :: p

  integer :: i, j, ifail, ielm_out
  real*8 :: x(3), R, Z, R_out, Z_out, Phi, tol, s, s_out, t_out
  real*8 :: m, e !< mean, rms error
  character*8 :: n_s, tol_s

  s = 0.d0
  if (present(smoothing)) s = smoothing
  call prepare_mumps_par(node_list, element_list, p, smoothing=s)

  do j=1,size(n)
    write(n_s, '(i8)') n(j)

    allocate(particles(n(j)))
    call find_RZ(node_list,element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)
    call initialise_particles(particles, node_list, element_list, rng)
    particles(:)%weight = volume/real(n(j))

    call project_particles(node_list, element_list, p, particles, 1)

    ! calculate tolerance based on number of points (to verify scaling)
    tol = 0.d0
    write(tol_s, '(g8.2)') tol
    ! test rms
    call elements_mean_rms(node_list, element_list, f_1, m, e)
    call assert_equals(1.d0, m, tol, 'mean value 1 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']')
    call assert_equals(0.d0, e, tol, 'rms value 0 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']')
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

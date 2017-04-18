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
use projection_spec, only: elements_mean_rms, f_1
use fruit
implicit none
include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure

logical, parameter :: write_proj_output = .true.

contains

subroutine test_square_10_10_pcg
  call project_n_square_10_10(10, pcg32_rng())
  call project_n_square_10_10(100, pcg32_rng())
  call project_n_square_10_10(1000, pcg32_rng())
  call project_n_square_10_10(10000, pcg32_rng())
  call project_n_square_10_10(100000, pcg32_rng())
end subroutine test_square_10_10_pcg

subroutine test_square_10_10_sob
  call project_n_square_10_10(10, sobseq_rng())
  call project_n_square_10_10(100, sobseq_rng())
  call project_n_square_10_10(1000, sobseq_rng())
  call project_n_square_10_10(10000, sobseq_rng())
  call project_n_square_10_10(100000, sobseq_rng())
end subroutine test_square_10_10_sob



subroutine project_n_square_10_10(n, rng)
  use mod_rng
  use constants, only: TWOPI
  use mod_sampling
  integer, intent(in) :: n
  class(type_rng), intent(in) :: rng
  integer, parameter :: n_R = 10, n_Z = 10
  type(particle_kinetic), dimension(:), allocatable :: particles
  type(DMUMPS_STRUC) :: p
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer :: i, ifail
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  real*8 :: x(3), R, Z, R_out, Z_out, Phi, tol
  real*8 :: m, e !< mean, rms error
  character*8 :: n_s, tol_s
  class(type_rng), allocatable :: my_rng

  allocate(my_rng, source=rng)
  call my_rng%initialize(n_dims=3, seed=1231789264, n_streams=1, i_stream=1, ierr=ifail)

  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)

  write(n_s, '(i8)') n

  allocate(particles(n))
  do i=1,size(particles,1)
    call my_rng%next(x)
    particles(i)%weight = TWOPI/real(size(particles,1))
    call transform_uniform_cylindrical(x, [R_geo-amin,R_geo+amin], [Z_geo-amin,Z_geo+amin],&
        [0.d0, TWOPI], R, Z, Phi)
    call find_RZ(node_list,element_list,R,Z,R_out,Z_out,&
        particles(i)%i_elm,particles(i)%st(1),particles(i)%st(2),ifail)
  end do

  call prepare_mumps_par(node_list, element_list, p, smoothing=0.d0)
  call project_particles(node_list, element_list, p, particles, 1)
  p%JOB=-2
  call DMUMPS(p)

  ! calculate tolerance based on number of points (to verify scaling)
  tol = 0.d0
  write(tol_s, '(g8.2)') tol
  ! test rms
  call elements_mean_rms(node_list, element_list, f_1, m, e)
  call assert_equals(0.d0, m, tol, 'rms value 0 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']')
  call assert_equals(1.d0, e, tol, 'mean value 1 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']')
  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, &
      filename='test_'//trim(adjustl(n_s))//'_square_10_10.h5', n_fields=1, time=0.d0)
  end if
end subroutine project_n_square_10_10
end module particle_projection_spec

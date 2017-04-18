!> This module contains some testcases for projecting particles, ensuring
!> that the projection matrix, RHS and MUMPS work for these cases.
!>
!> It contains tests of projecting zero, 1, x, xy onto a square or circular grid
module projection_spec
use mod_project_particles
use data_structure
use mod_particle_types
use fruit
implicit none
include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure

logical, parameter :: write_proj_output = .true.

contains

!> Actions to perform before any of these tests
subroutine setup_projection_spec
  call initialise_basis
end subroutine setup_projection_spec

!> Functions to project
function f_0(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_0
  f_0 = 0.d0
end function f_0
function f_1(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_1
  f_1 = 1.d0
end function f_1
function f_R(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_R
  f_R = R
end function f_R
function f_RZ(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_RZ
  f_RZ = R*Z
end function f_RZ
function f_R4(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_R4
  f_R4 = R**4
end function f_R4


subroutine default_square_grid(node_list, element_list, n)
  type(type_node_list), intent(out) :: node_list
  type(type_element_list), intent(out) :: element_list
  integer, intent(in) :: n !< number of nodes in each dimension
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n, n, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
end subroutine


!> Project zeros onto a square grid
subroutine test_project_0_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_0, 0.d0, 0.d0, '0_square_10_10')
end subroutine test_project_0_square_10_10

!> Project one onto a square grid
subroutine test_project_1_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer :: index(1:10*10*4), i, j

  call default_square_grid(node_list, element_list, 10)

  ! include a grid_bezier_square test here
  ! verify that no node shares the same index and all indices are used exactly once
  index = 0
  do i=1,node_list%n_nodes
    do j=1,n_order+1
      index(node_list%node(i)%index(j)) = index(node_list%node(i)%index(j)) + 1
    end do
  end do
  call assert_equals(size(index,1), count(index .gt. 0), 'all indices must be used')
  call assert_equals(0, count(index .gt. 1), 'no duplicate indices in this grid')
  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_square_10_10')
end subroutine test_project_1_square_10_10

subroutine test_project_1_polar_101_32
  use phys_module
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call initialise_basis
  call preset_parameters()
  fbnd(1) = 2.d0
  fbnd(2:4) = 0.d0
  mf = 0
  n_radial = 101
  n_pol = 32
  R_geo = 1.5
  Z_geo = 0.0
  amin = 1.0

  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
    node_list, element_list)

  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_polar_101_32')
end subroutine test_project_1_polar_101_32


!> Project R onto a square grid
subroutine test_project_R_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_R, 26.d0/24.d0, 0.d0, 'R_square_10_10')
end subroutine test_project_R_square_10_10

!> Project RZ onto a square grid
subroutine test_project_RZ_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_RZ, 0.d0, 0.d0, 'RZ_square_10_10')
end subroutine test_project_RZ_square_10_10

subroutine test_project_R4_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_R4, 1.89583333333333, 0.d0, 'R4_square_10_10', rms_tol=3d-6)
end subroutine test_project_R4_square_10_10

subroutine test_project_R4_square_n
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_max = 40
  integer :: i
  character(len=13) :: s
  real*8 :: tol
  do i=3,n_max,4
    write(s,'(A,i0.3)') "R4_square_", i
    call default_square_grid(node_list, element_list, i)
    tol = 0.06d0/(real(i)**4)
    call project_f_with_assert_and_write(node_list, element_list, f_R4, 1.89583333333333d0, 0.d0, s, rms_tol=tol)
  end do
end subroutine test_project_R4_square_n




!> Project a function f onto grid in node_list and element_list
!> and test for mean and RMS value. Optionally write to file for visual inspection.
subroutine project_f_with_assert_and_write(node_list, element_list, f, mean, RMS, name, rms_tol)
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, external :: f
  real*8, intent(in) :: mean
  real*8, intent(in) :: rms
  real*8, intent(in), optional :: rms_tol
  character(len=*), intent(in) :: name
  real*8 :: m, e, my_rms_tol
  call project_f(node_list, element_list, f)
  ! test rms
  my_rms_tol = 1d-12
  if (present(rms_tol)) my_rms_tol = rms_tol
  call elements_mean_rms(node_list, element_list, f, m, e)
  call assert_equals(mean, m, 1d-12, 'mean value M')
  call assert_equals(RMS, e, my_rms_tol, 'rms error ok')
  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, &
      filename=name//'.h5', n_fields=1, time=0.d0)
  end if
end subroutine project_f_with_assert_and_write



subroutine test_projection_matrix_square_2_2
  use basis_at_gaussian
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  type(DMUMPS_STRUC) :: p
  integer :: i, j

  integer, parameter :: n_R = 2, n_Z = 2
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0

  real*8, parameter :: ref(16) = [0.100816,0.0159184,0.0142177,0.0022449,& ! index 1
      0.0477551,-0.0110544,0.00673469,-0.00155896,& ! index 2
      0.034898,0.0055102,-0.00840136,-0.00132653,& ! index 3 (but node 4, because index is switched with matrix order)
      0.0165306,-0.00382653,-0.00397959,0.000921202] ! index 4 (but node 3)
  real*8, parameter :: tol = 1d-6
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  call prepare_mumps_par(node_list, element_list, p, smoothing=0d0, skip_factorisation=.true.)
  do i=1,size(p%irn)
    if (p%irn(i) == 1) call assert_equals(ref(p%jcn(i)), p%A(i), tol, 'matrix element must match reference')
  end do
end subroutine test_projection_matrix_square_2_2




! ****************
! Helper functions
! ****************

!> Project a function onto the JOREK elements
subroutine project_f(node_list, element_list, f)
  use basis_at_gaussian
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, external :: f
  type(DMUMPS_STRUC) :: p
  integer :: i, j, k, m, index, i_elm, inode, ms, mt
  real*8, dimension(n_gauss,n_gauss) :: x_g, y_g, x_s, x_t, y_s, y_t
  real*8 :: wst, xjac, v
  type(type_node) :: nodes(4)
  type(type_element) :: element

  call prepare_mumps_par(node_list, element_list, p, smoothing=0d0)

  ! Project manually
  p%JOB = 3
  p%icntl(21) = 0 ! solution is available only on host
  p%icntl(4)  = 1 ! print only errors
  ! Setup RHS by integrating manually
  p%rhs = 0.d0

  do i_elm=1,element_list%n_elements
    element = element_list%element(i_elm)
    do m=1,n_vertex_max
      nodes(m) = node_list%node(element%vertex(m))
    enddo

    ! Set up gauss points in this element
    x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
    do i=1,n_vertex_max
      do j=1,n_order+1
        do ms=1, n_gauss
          do mt=1, n_gauss
            x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
            y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

            x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
            x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

            y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
            y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)
          enddo
        enddo
      enddo
    enddo

    ! Perform gauss integration of RHS
    do ms=1, n_gauss
      do mt=1, n_gauss
        wst = wgauss(ms)*wgauss(mt)
        xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

        do i=1,n_vertex_max
          do j=1,n_order+1
            index = nodes(i)%index(j)

            v   = h(i,j,ms,mt)  * element%size(i,j)
            p%rhs(index) = p%rhs(index) + f(x_g(ms,mt), y_g(ms,mt)) * v * xjac * x_g(ms,mt) * wst
          enddo
        enddo
      enddo
    enddo
  enddo
  call DMUMPS(p)

  do i=1,node_list%n_nodes
    do k=1,n_order+1
      index = node_list%node(i)%index(k)
      node_list%node(i)%values(1,k,1) = p%rhs(index)
      write(*,*) index, node_list%node(i)%values(1,k,1), p%rhs(index)
    enddo
  enddo

  p%JOB=-2
  call DMUMPS(p)
end subroutine project_f



subroutine elements_mean_rms(node_list, element_list, f, mean, rms)
  use basis_at_gaussian
  use constants, only: TWOPI
  use mod_interp4
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8, external :: f
  real*8, intent(out) :: mean, rms

  integer :: i_elm, m, i, j, ms, mt
  type(type_element) :: element
  type(type_node) :: nodes(4)
  real*8, dimension(n_gauss,n_gauss) :: x_g, x_s, x_t, y_g, y_s, y_t
  real*8 :: my_ref, wst, volume, xjac
  real*8, dimension(1) :: P, P_s, P_t, P_st, P_ss, P_tt

  call initialise_basis

  volume = 0.d0
  mean = 0.d0
  rms = 0.d0

  do i_elm=1,element_list%n_elements

    element = element_list%element(i_elm)
    do m=1,n_vertex_max
      nodes(m) = node_list%node(element%vertex(m))
    enddo

    ! Set up gauss points in this element
    x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
    do i=1,n_vertex_max
      do j=1,n_order+1
        do ms=1, n_gauss
          do mt=1, n_gauss
            x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
            y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

            x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
            x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)

            y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
            y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)
          enddo
        enddo
      enddo
    enddo

    ! Perform gauss integration of LHS
    do ms=1, n_gauss
      do mt=1, n_gauss
        wst = wgauss(ms)*wgauss(mt)
        xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        volume = volume + TWOPI * x_g(ms,mt) * xjac * wst

        ! calculate contribution to integral of this point
        call interp(node_list, element_list, i_elm, 1, 1, Xgauss(ms), Xgauss(mt), P, P_s, P_t, P_st, P_ss, P_tt)
        rms = rms + (P(1)-f(x_g(ms,mt),y_g(ms,mt)))**2 * xjac * TWOPI * x_g(ms,mt) * wst
        mean = mean + P(1) * xjac * TWOPI * x_g(ms,mt) * wst
      enddo
    enddo
  enddo

  rms = sqrt(rms / volume)
  mean = mean/volume
end subroutine elements_mean_rms
end module projection_spec

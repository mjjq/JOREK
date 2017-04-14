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

!> Project zeros onto a square grid
subroutine test_project_0_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_R = 10, n_Z = 10
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  call project_0(node_list, element_list, '0_square_10_10')
end subroutine test_project_0_square_10_10
!> Test projecting a constant function (0) onto some grid
!> Calculate the rhs by integrating this function.
!> The result should be a constant 0 everywhere. Verify this.
subroutine project_0(node_list, element_list, name)
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*), intent(in) :: name
  integer :: i
  call project_f(node_list, element_list, f_0)
  ! test a few positions
  do i=1,20
    call assert_equals(0.d0, node_list%node(nint(real(i)/20.d0*real(node_list%n_nodes)))%values(1,1,1), 1d-9, 'value must be 0')
  enddo
  ! test rms
  call assert_equals(0.d0, elements_rms(node_list, element_list, 0.d0), 1d-12, 'rms value 0')
  call assert_equals(0.d0, elements_mean(node_list, element_list), 1d-12, 'mean value 0')
  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, &
      filename=name//'.h5', n_fields=1, time=0.d0)
  end if
end subroutine project_0
function f_0(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_0
  f_0 = 0
end function f_0




!> Project one onto a square grid
subroutine test_project_1_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_R = 10, n_Z = 10
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  integer :: index(1:10*10*4), i, j
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)

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
  call project_1(node_list, element_list, '1_square_10_10')
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

  call project_1(node_list, element_list, '1_polar_101_32')
end subroutine test_project_1_polar_101_32


!> Test projecting a constant function (1) onto some grid
!> Calculate the rhs by integrating this function.
!> The result should be a constant 1 everywhere. Verify this.
subroutine project_1(node_list, element_list, name)
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*), intent(in) :: name
  integer :: i
  call project_f(node_list, element_list, f_1)
  ! test a few positions
  do i=1,20
    call assert_equals(1.d0, node_list%node(nint(real(i)/20.d0*real(node_list%n_nodes)))%values(1,1,1), 1d-9, 'value must be 1')
  enddo
  ! test rms
  call assert_equals(0.d0, elements_rms(node_list, element_list, 1.d0), 7d-5, 'rms value 0')
  call assert_equals(1.d0, elements_mean(node_list, element_list), 1d-6, 'mean value 1')
  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, &
      filename=name//'.h5', n_fields=1, time=0.d0)
  end if
end subroutine project_1
function f_1(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_1
  f_1 = 1
end function f_1


!> Test projecting r onto some grid
!> Calculate the rhs by integrating this function.
!> check the mean and rms value (calculated here:)
!> mean = int r^2 dr dz dphi / volume = 26/24 = M
!> ms = int (r-M)^2 r dr = 11/144
!> rms = sqrt(ms) = sqrt(11/144) \approx 0.276385
!> WARNING: valid for specific square grid only
subroutine project_R(node_list, element_list, name)
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*), intent(in) :: name
  integer :: i
  ! WARNING: valid for specific square grid only
  real*8, parameter :: M = 26.d0/24.d0
  real*8, parameter :: RMS = sqrt(11.d0/144.d0)
  call project_f(node_list, element_list, f_r)
  ! test rms
  call assert_equals(M, elements_mean(node_list, element_list), 1d-12, 'mean value 26/24')
  call assert_equals(RMS, elements_rms(node_list, element_list, M), 1d-12, 'rms value ok')
  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, &
      filename=name//'.h5', n_fields=1, time=0.d0)
  end if
end subroutine project_R
function f_R(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_R
  f_R = R
end function f_R
subroutine test_project_R_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_R = 10, n_Z = 10
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  call project_R(node_list, element_list, 'R_square_10_10')
end subroutine test_project_R_square_10_10
subroutine test_project_R_square_50_50
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_R = 50, n_Z = 50
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  call project_R(node_list, element_list, 'R_square_50_50')
end subroutine test_project_R_square_50_50


!> Test projecting rz onto some grid
!> Calculate the rhs by integrating this function.
!> check the mean and rms value (calculated here:)
!> mean = int z r^2 dr dz dphi / volume = 0 (because odd)
!> ms = int (r-M)^2 r dr = 5/48
!> rms = sqrt(ms) = sqrt(5/48) \approx 0.322749
!> WARNING: valid for specific square grid only
subroutine project_RZ(node_list, element_list, name)
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*), intent(in) :: name
  integer :: i
  ! WARNING: valid for specific square grid only
  real*8, parameter :: M = 0.d0
  real*8, parameter :: RMS = sqrt(5.d0/48.d0)
  call project_f(node_list, element_list, f_rz)
  ! test rms
  call assert_equals(M, elements_mean(node_list, element_list), 1d-12, 'mean value 0')
  call assert_equals(RMS, elements_rms(node_list, element_list, M), 1d-12, 'rms value ok')
  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, &
      filename=name//'.h5', n_fields=1, time=0.d0)
  end if
end subroutine project_RZ
function f_RZ(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_RZ
  f_RZ = R*Z
end function f_RZ
subroutine test_project_RZ_square_10_10
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_R = 10, n_Z = 10
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  call project_RZ(node_list, element_list, 'RZ_square_10_10')
end subroutine test_project_RZ_square_10_10
subroutine test_project_RZ_square_50_50
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_R = 50, n_Z = 50
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  call initialise_basis
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  call project_RZ(node_list, element_list, 'RZ_square_50_50')
end subroutine test_project_RZ_square_50_50



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
    if (p%irn(i) == 1) write(*,*) p%jcn(i), p%A(i)
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

  call prepare_mumps_par(node_list, element_list, p, smoothing=0d-3)

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
  !write(*,*) p%rhs
  write(*,*) p%rhs(1:30)

  call DMUMPS(p)
  write(*,*) "solved rhs"
  write(*,*) p%rhs(1:30)

  do i=1,node_list%n_nodes
    do k=1,n_order+1
      index = node_list%node(i)%index(k)
      node_list%node(i)%values(1,k,1) = p%rhs(index)
    enddo
  enddo

  p%JOB=-2
  call DMUMPS(p)
end subroutine project_f



!> For ease of calling. Does the work twice if you need both mean and rms though..
function elements_mean(node_list, element_list)
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8 :: elements_mean
  real*8 :: tmp
  call elements_mean_rms(node_list, element_list, elements_mean, tmp)
end function elements_mean
!> For ease of calling. Does the work twice if you need both mean and rms though..
function elements_rms(node_list, element_list, ref)
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8, intent(in), optional :: ref
  real*8 :: elements_rms
  real*8 :: tmp
  if (present(ref)) then
    call elements_mean_rms(node_list, element_list, tmp, elements_rms, ref)
  else
    call elements_mean_rms(node_list, element_list, tmp, elements_rms)
  end if
end function elements_rms

subroutine elements_mean_rms(node_list, element_list, mean, rms, ref)
  use basis_at_gaussian
  use constants, only: TWOPI
  use mod_interp4
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8, intent(in), optional :: ref
  real*8, intent(out) :: mean, rms

  integer :: i_elm, m, i, j, ms, mt
  type(type_element) :: element
  type(type_node) :: nodes(4)
  real*8, dimension(n_gauss,n_gauss) :: x_g, x_s, x_t, y_g, y_s, y_t
  real*8 :: my_ref, wst, volume, xjac, P(1)

  my_ref = 0.d0
  if (present(ref)) my_ref = ref
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
        call interp4(node_list, element_list, i_elm, [1], 1, Xgauss(ms), Xgauss(mt), 0.d0, P)
        rms = rms + (P(1)-my_ref)**2 * xjac * TWOPI * x_g(ms,mt) * wst
        mean = mean + P(1) * xjac * TWOPI * x_g(ms,mt) * wst
      enddo
    enddo
  enddo

  rms = sqrt(rms / volume)
  mean = mean/volume
end subroutine elements_mean_rms
end module projection_spec

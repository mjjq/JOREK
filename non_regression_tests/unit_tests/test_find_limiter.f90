!> Verify the generation of a few simple grids
module test_find_limiter
use fruit
use data_structure
use mod_boundary
implicit none
contains

!> Test if we can find the limiter on a square grid with a linear variation of psi.
!> the grid consists only of a single element.
subroutine test_square_grid
  type(type_node_list) :: node_list
  type(type_element_list) :: element_list
  type (type_bnd_node_list) :: bnd_node_list
  type (type_bnd_element_list) :: bnd_elm_list
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  integer :: i, j
  real*8 :: psi_lim, R_lim, Z_lim
  logical, parameter :: write_psi = .false.
  
  ! for find_RZ and interp
  real*8 :: R_out, Z_out, s, t
  integer :: i_elm, ifail
  real*8, dimension(1) :: P, P_s, P_t, P_st, P_tt, P_ss

  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(2, 2, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)

  ! Fill in some analytical values, only for the values and derivatives
  ! the cross-derivatives are not taken into account by find_limiter and psi_minmax

  do i=1,node_list%n_nodes
    node_list%node(i)%values(1,1,1) = 0.d0 ! value at node
    node_list%node(i)%values(1,3,1) = 0.d0 ! derivative in 'z'-direction
    node_list%node(i)%values(1,4,1) = 0.d0 ! skip this one, it's the cross derivative and not used

    if (node_list%node(i)%x(1,2) .le. 0.d0) then
      node_list%node(i)%values(1,2,1) = -1.d0 ! derivative in 'r'-direction (over element size)
    else
      node_list%node(i)%values(1,2,1) = -0.5d0 ! derivative in 'r'-direction (over element size)
    end if
  end do

  ! create boundary
  call boundary_from_grid(node_list,element_list,bnd_node_list,bnd_elm_list,infos=.false.)

  call find_limiter(0, node_list, element_list, bnd_elm_list, psi_lim, R_lim, Z_lim)

  ! now we expect to have a minimum at RZ = (1.0, -0.5)
  call assert_equals(1.d0, R_lim, 1d-12, 'R_limiter')
  call assert_equals(-0.5d0, Z_lim, 1d-12, 'Z_limiter')
  call assert_equals(-0.5d0, psi_lim, 1d-12, 'psi_limiter')

  call find_RZ(node_list, element_list,R_lim,Z_lim,R_out,Z_out,i_elm,s,t,ifail)
  call interp(node_list, element_list, i_elm, 1, 1, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
  call assert_equals(psi_lim, P(1), 1d-12, 'psi_limiter equal to psi at R_lim, Z_lim')

  if (write_psi) then
    open(unit=12,file='psi.txt')
    do i=1,100
      do j=1,100
        s = real(i-1,8)/99.d0
        t = real(j-1,8)/99.d0
        call interp(node_list, element_list, i_elm, 1, 1, s,t, P, P_s, P_t, P_st, P_ss, P_tt)
        write(12,*) s,t,P
      end do
    end do
    close(12)
  end if
end subroutine test_square_grid
end module test_find_limiter

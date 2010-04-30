subroutine boundary_check
use data_structure
use phys_module
use vacuum_response_module
use nodes_elements

implicit none

integer :: i_harm, ibnd, i_elm, iv1, iv2, inode1, inode2, inode1a, inode2a, kbnd, lbnd
integer :: index_basis_bnd, index_basis2_bnd, n_points, i, j, ij_node1, ij_node2, ij_node, kl_node, ip
real*8  :: s_out, t_out, B_n, B_tan, B_tan2, B_tan_v, xjac, grad_ss, grad_st, B_n_check
real*8  :: H1(2,2), H1_s(2,2), H1_ss(2,2)
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt

write(*,*) '************************************'
write(*,*) '*    check boundary conditions     *'
write(*,*) '************************************'

i_harm = 2

!  type type_boundary                                  ! type definition for one boundary element (1D element)
!    integer :: vertex(2)                              ! the nodes of the corners
!    integer :: direction(2,2)                         ! indicates which direction of the nodes is along the boundary (2 or 3)
!    integer :: element                                ! boundary element is part of this element
!    integer :: side                                   ! boundary element corresponds to this side of the originating element
!    real*8  :: size(2,2)                              ! the size of the vectors at each vertex of the element : size(vertex,order)
!  endtype type_boundary

n_points = 11

if (freeboundary) then

  do ibnd = 1, boundary_list%n_boundary

    do ip = 1, n_points

      i_elm = boundary_list%boundary(ibnd)%element
      iv1   = boundary_list%boundary(ibnd)%side
      iv2   = mod(iv1,4)+1

      inode1 = element_list%element(i_elm)%vertex(iv1)
      inode2 = element_list%element(i_elm)%vertex(iv2)

      inode1a = boundary_list%boundary(ibnd)%vertex(1)
      inode2a = boundary_list%boundary(ibnd)%vertex(2)

      if ((iv1 .eq. 1) .and. (iv2 .eq. 2)) then
        t_out = 0.d0
        s_out = float(ip-1)/float(n_points)
        call basisfunctions1(s_out,H1,H1_s,H1_ss)
      elseif ((iv1 .eq. 2) .and. (iv2 .eq. 3)) then
        t_out = float(ip-1)/float(n_points)
        s_out = 1.d0
        call basisfunctions1(t_out,H1,H1_s,H1_ss)
      elseif ((iv1 .eq. 3) .and. (iv2 .eq. 4)) then
        t_out = 1.d0
        s_out = float(ip-1)/float(n_points)
        call basisfunctions1(s_out,H1,H1_s,H1_ss)
      elseif ((iv1 .eq. 4) .and. (iv2 .eq. 1)) then
        t_out = float(ip-1)/float(n_points)
        s_out = 0.d0
        call basisfunctions1(t_out,H1,H1_s,H1_ss)
      endif

      call interp_RZ(node_list,element_list,i_elm,s_out,t_out,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

      call interp(node_list,element_list,i_elm,1,i_harm,s_out,t_out,P,P_s,P_t,P_st,P_ss,P_tt)

      xjac = R_s * Z_t - R_t * Z_s

      grad_ss =   (R_t**2 + Z_t**2)       / xjac**2
      grad_st = - (R_s * R_t + Z_s * Z_t) / xjac**2

      B_tan = (P_s * grad_ss + P_t * grad_st) / (R**2 * sqrt(grad_ss))
      B_n   = P_t / (R * xjac * sqrt(grad_ss))

!--------------------------------------------------------------check whether information from boundary gives same result as element-wise
      ij_node1   = boundary_list%boundary(ibnd)%vertex(1)
      ij_node2   = boundary_list%boundary(ibnd)%vertex(2)

      B_n_check =             node_list%node(ij_node1)%values(i_harm,1,1) * boundary_list%boundary(ibnd)%size(1,1) * H1_s(1,1)
      B_n_check = B_n_check + node_list%node(ij_node1)%values(i_harm,3,1) * boundary_list%boundary(ibnd)%size(1,2) * H1_s(1,2)
      B_n_check = B_n_check + node_list%node(ij_node2)%values(i_harm,1,1) * boundary_list%boundary(ibnd)%size(2,1) * H1_s(2,1)
      B_n_check = B_n_check + node_list%node(ij_node2)%values(i_harm,3,1) * boundary_list%boundary(ibnd)%size(2,2) * H1_s(2,2)

      B_n_check = B_n_check / (R * xjac * sqrt(grad_ss))

      B_tan_v = 0.d0

      do i = 1, 2                                                                  ! two vertices of a boundary element

        ij_node = boundary_list%boundary(ibnd)%vertex(i)

        do j = 1, 2                                                                ! two basisfunctions at each vertex

          index_basis_bnd = 2*mod(ibnd+i-2,boundary_list%n_boundary) + (j-1) + 1   ! the index in the vacuum_response matrix (DANGEROUS ASSUMPTION!)

          do kbnd=1,boundary_list%n_boundary                                       ! sum over all boudary nodes

            kl_node = boundary_list%boundary(kbnd)%vertex(1)                       ! dangerous : can give same node twice if orientation changes

            do lbnd = 1, 2

              index_basis2_bnd = 2*(kbnd-1) +   (lbnd-1) + 1                       ! the index in the vacuum_response matrix

              P_t = node_list%node(kl_node)%values(i_harm,2*lbnd-1,1)              ! to be changed now only using index 1 and 3

!              B_tan_v = B_tan_v + vacuum_response(index_basis_bnd,index_basis2_bnd,i_harm) * P_t * boundary_list%boundary(ibnd)%size(i,j) * H1_s(i,j)

              B_tan_v = B_tan_v + vacuum_response(index_basis2_bnd,index_basis_bnd,i_harm) * P_t * boundary_list%boundary(ibnd)%size(i,j) * H1_s(i,j)

            enddo
          enddo

        enddo
      enddo

      write(*,'(i5,8e14.6)') n_points*(ibnd-1)+ip,B_n,B_tan,B_tan_v

    enddo

  enddo

endif

return
end

subroutine update_neighbours(element_list,node_list)
use data_structure
use mod_neighbours
use mod_element_rtree
use mod_qsort
implicit none

type (type_element_list) :: element_list
type (type_node_list)    :: node_list
type (type_element)      :: elm_i, elm_j
integer                  :: inb_i, inb_j, i, j, k
integer                  :: i_elm, j_elm, i_node1, i_node2
real*8                   :: s_i, t_i, R_i, Rs_i, Rt_i, Rst_i, Rss_i, Rtt_i, Z_i, Zs_i, Zt_i, Zst_i,Zss_i,Ztt_i
real*8                   :: s_j, t_j, R_j, Rs_j, Rt_j, Rst_j, Rss_j, Rtt_j, Z_j, Zs_j, Zt_j, Zst_j,Zss_j,Ztt_j
integer, dimension(:), allocatable :: i_nearby

call populate_element_rtree(node_list, element_list)

!$omp parallel do default(private) &
!$omp   shared(element_list,node_list)
do i=1, element_list%n_elements
  call nearby_elements(node_list, element_list, i, i_nearby)
  do k=1,size(i_nearby,1)
    j = i_nearby(k)
    if (i .eq. j) cycle

    if  (neighbours(node_list, element_list%element(i), element_list%element(j),inb_i,inb_j)) then

       !write(*,*) 'found neighbours : ',i,j,inb_i,inb_j

        element_list%element(i)%neighbours(inb_i) = j

      if (abs(inb_i-inb_j) .eq. 2) then

        if (inb_i .eq. 1) then                     ! must check for clockwise or anti-clockwise    ! should also take into account scale lengts between old and new (s,t)
          element_list%element(i)%transform(inb_i,1,:) = (/0., 1., 0. /)                ! s_new = s_old
          element_list%element(i)%transform(inb_i,2,:) = (/1., 0., 1. /)                ! t_new = 1. + t_old
          s_i = 0.371
          t_i = 0.
        elseif (inb_i .eq. 3) then
          element_list%element(i)%transform(inb_i,1,:) = (/ 0., 1., 0. /)               ! s_new = s_old
          element_list%element(i)%transform(inb_i,2,:) = (/-1., 0., 1. /)               ! t_new = t_old - 1.
          s_i = 0.371
          t_i = 1.
        elseif (inb_i .eq. 2) then
          element_list%element(i)%transform(inb_i,1,:) = (/-1., 1., 0. /)               ! s_new = s_old - 1.
          element_list%element(i)%transform(inb_i,2,:) = (/ 0., 0., 1. /)                ! t_new = t_old
          s_i = 1.
          t_i = 0.371
        elseif (inb_i .eq. 4) then
          element_list%element(i)%transform(inb_i,1,:) = (/ 1., 1., 0. /)                ! s_new = 1. + sold
          element_list%element(i)%transform(inb_i,2,:) = (/ 0., 0., 1. /)                ! t_new = t_old
          s_i = 0
          t_i = 0.371
        endif

      i_elm = i
      j_elm = element_list%element(i)%neighbours(inb_i)

      s_j = element_list%element(i)%transform(inb_i,1,1) +  element_list%element(i)%transform(inb_i,1,2) * s_i +   element_list%element(i)%transform(inb_i,1,3)* t_i
      t_j = element_list%element(i)%transform(inb_i,2,1) +  element_list%element(i)%transform(inb_i,2,2) * s_i +   element_list%element(i)%transform(inb_i,2,3)* t_i

      call interp_RZ(node_list,element_list,i_elm,s_i,t_i,R_i,Rs_i,Rt_i,Rst_i,Rss_i,Rtt_i,Z_i,Zs_i,Zt_i,Zst_i,Zss_i,Ztt_i)
      call interp_RZ(node_list,element_list,j_elm,s_j,t_j,R_j,Rs_j,Rt_j,Rst_j,Rss_j,Rtt_j,Z_j,Zs_j,Zt_j,Zst_j,Zss_j,Ztt_j)

!      write(*,'(A,2i5,8f9.4)') 'i_elm ',i_elm,inb_i,s_i,t_i,R_i,Z_i
!      write(*,'(A,2i5,8f9.4)') 'j_elm ',j_elm,inb_j,s_j,t_j,R_j,Z_j

      if ( abs(R_i-R_j)+ abs(Z_i-Z_j) .gt. 1.e-8) then
        write(*,'(A,2i5,8f16.10)') 'PROBLEM : ',i_elm,j_elm,R_i,R_j,Z_i,Z_j
      endif

      else
        !write(*,*) 'not yet implemented '
      endif

    endif
  enddo
enddo
!$omp end parallel do

do i=1, element_list%n_elements

  do  j= 1, n_vertex_max

    i_node1 = element_list%element(i)%vertex(j)
    i_node2 = element_list%element(i)%vertex(mod(j,4)+1)

    if (norm2(node_list%node(i_node1)%x(1,1:2) - node_list%node(i_node2)%x(1,1:2)) .lt. 1d-8) then
      element_list%element(i)%neighbours(j) = -1
      do k=1,n_vertex_max
         if(element_list%element(i)%neighbours(k).eq.0) then
            element_list%element(i)%neighbours(k) = -1
         end if
      end do
    endif

  enddo

enddo


!do i=1, element_list%n_elements
!  write(*,'(12i6)') i, element_list%element(i)%vertex, element_list%element(i)%neighbours
!  do  j= 1, 4
!    write(*,'(i5,6f12.4)') j,element_list%element(i)%transform(j,1,:),element_list%element(i)%transform(j,2,:)
!  enddo
!enddo
end

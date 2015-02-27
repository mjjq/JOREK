subroutine update_neighbours(element_list,node_list)
use data_structure
implicit none

type (type_element_list) :: element_list
type (type_node_list)    :: node_list
type (type_element)      :: elm_i, elm_j
integer                  :: inb_i, inb_j, i, j
integer                  :: i_elm, j_elm
real*8                   :: s_i, t_i, R_i, Rs_i, Rt_i, Z_i, Zs_i, Zt_i, s_j, t_j, R_j, Rs_j, Rt_j, Z_j, Zs_j, Zt_j
logical, external        :: neighbours

write(8,*) 'updating neighbours'

do i=1, element_list%n_elements

  do j=1, element_list%n_elements

    if (i .ne. j) then

    if  (neighbours( element_list%element(i), element_list%element(j),inb_i,inb_j)) then

       !write(*,*) 'found neighbours : ',i,j,inb_i,inb_j

        element_list%element(i)%neighbours(inb_i) = j

      if (abs(inb_i-inb_j) .eq. 2) then

        if (inb_i .eq. 1) then                     ! must check for clockwise or anti-clockwise
          element_list%element(i)%transform(inb_i,1,:) = (/0., 1., 0. /)
          element_list%element(i)%transform(inb_i,2,:) = (/1., 0., 0. /)
          s_i = 0.371
          t_i = 0.
        elseif (inb_i .eq. 3) then
          element_list%element(i)%transform(inb_i,1,:) = (/0., 1., 0. /)
          element_list%element(i)%transform(inb_i,2,:) = (/0., 0., 0. /)
          s_i = 0.371
          t_i = 1.
        elseif (inb_i .eq. 2) then
          element_list%element(i)%transform(inb_i,1,:) = (/0., 0., 0. /)
          element_list%element(i)%transform(inb_i,2,:) = (/0., 0., 1. /)
          s_i = 1.
          t_i = 0.371
        elseif (inb_i .eq. 4) then
          element_list%element(i)%transform(inb_i,1,:) = (/1., 0., 0. /)
          element_list%element(i)%transform(inb_i,2,:) = (/0., 0., 1. /)
          s_i = 0
          t_i = 0.371
        endif

      else
        write(*,*) 'not yet implemented'
      endif

      i_elm = i
      j_elm = element_list%element(i)%neighbours(inb_i)

      s_j = element_list%element(i)%transform(inb_i,1,1) +  element_list%element(i)%transform(inb_i,1,2) * s_i +   element_list%element(i)%transform(inb_i,1,3)* t_i
      t_j = element_list%element(i)%transform(inb_i,2,1) +  element_list%element(i)%transform(inb_i,2,2) * s_i +   element_list%element(i)%transform(inb_i,2,3)* t_i

      call interp_RZ(node_list,element_list,i_elm,s_i,t_i,R_i,Rs_i,Rt_i,Z_i,Zs_i,Zt_i)
      call interp_RZ(node_list,element_list,j_elm,s_j,t_j,R_j,Rs_j,Rt_j,Z_j,Zs_j,Zt_j)

!      write(*,'(A,2i5,8f9.4)') 'i_elm ',i_elm,inb_i,s_i,t_i,R_i,Z_i
!      write(*,'(A,2i5,8f9.4)') 'j_elm ',j_elm,inb_j,s_j,t_j,R_j,Z_j

      if ( abs(R_i-R_j)+ abs(Z_i-Z_j) .gt. 1.e-8) then
        write(*,'(A,2i5,8f16.10)') 'PROBLEM : ',i_elm,j_elm,R_i,R_j,Z_i,Z_j
      endif

    endif
    endif

  enddo
enddo

return
end
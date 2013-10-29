subroutine update_values(my_id,node_list,RHS)
!-----------------------------------------------------------------------
! subroutine adds the delta_values in RHS to the values in the node_list
!-----------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list) :: node_list
real*8  :: RHS(*)
integer :: my_id, i, j, k, in, index_node, index

if (my_id .eq. 0) then

  do i = 1, node_list%n_nodes

    do j=1,n_order+1

      index_node = node_list%node(i)%index(j)

      do k=1,n_var
        
!	if ((k .eq. 2) .or. (k .eq. 4)) then       !IF statement for simulations without flow
	
!	else                                       !IF statement for simulations without flow
	 
           do in=1,n_tor

             index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in

             if (index .gt. 0) then

               node_list%node(i)%values(in,j,k) = node_list%node(i)%values(in,j,k) + RHS(index)
               node_list%node(i)%deltas(in,j,k) = RHS(index)

             endif

           enddo
	       
!	endif                                      !IF statement for simulations without flow
	   
      enddo

    enddo

!   write(*,'(i5,20e12.4)') i,node_list%node(i)%values(1,:,2),node_list%node(i)%values(2,:,2)

  enddo

endif

call broadcast_nodes(my_id,node_list)

return
end

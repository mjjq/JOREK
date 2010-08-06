subroutine distribute_nodes_elements(my_id,n_cpu,node_list,element_list, &
                                    local_elms, n_local_elms, n_dof, index_min, index_max)
!---------------------------------------------------------------------------------------------
! subroutine divides the nodes (not their individual dof) over n_cpu equal parts
!            builds local_elms, contain all elements with at least one node with 
!            one index between index_min and index_max
!---------------------------------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: flux_list

integer :: local_elms(*)
integer :: my_id, n_cpu, n_dof, n_local_elms, index_total, inode
integer :: index_min(*), index_max(*), index_part, inext, i,j, k, iv,index1

logical :: elm_is_local
!integer, dimension(node_list%n_nodes) :: active_node
!integer                               :: n_active_nodes
if (my_id .eq.0) then
  write(*,*) '************************************'
  write(*,*) '*     distributing nodes           *'
  write(*,*) '************************************'
endif


 

  
 !call  Ref_Active_node( element_list,node_list ,active_node,n_active_nodes)
 !index_total = -1
 !do i=1,n_active_nodes
  !inode=active_node(i)
  !index_total = max(index_total,maxval(node_list%node(inode)%index))
! enddo

 index_total = -1
do inode=1,node_list%n_nodes
  index_total = max(index_total,maxval(node_list%node(inode)%index))
enddo
!stop
!write(*,*) ' n_elements  : ',my_id,element_list%n_elements
!write(*,*) ' n_nodes     : ',my_id,node_list%n_nodes
!write(*,*) ' index_total : ',my_id,index_total

index_min(1:n_cpu) = 0
index_max(1:n_cpu) = 0

!----------------------------- must really take into account the number of elements contributing to each node
index_min(1) = 1
do i=1,n_cpu
  index_max(i) = (i * index_total) / n_cpu
enddo
do i=2,n_cpu
  index_min(i) = index_max(i-1) + 1
enddo
if (my_id .eq. n_cpu-1) index_max(my_id+1) = index_total

!write(*,'(A,3i6)') ' index_min,index_max : ',my_id,index_min(my_id+1),index_max(my_id+1)
!write(*,'(A,3i6)') ' index_part          : ',my_id,index_part

n_dof           = index_total * n_tor * n_var

!----------------------------------------------- find the elements that have a local node
inext = 0

do i = 1, element_list%n_elements

  ELM_is_local = .false.

  do iv=1,n_vertex_max

    inode = element_list%element(i)%vertex(iv)

    do k=1, n_order+1

      if ( (node_list%node(inode)%index(k) .ge. index_min(my_id+1)) .and. &
           (node_list%node(inode)%index(k) .le. index_max(my_id+1)) ) then
        ELM_is_local = .true.
        goto 10
      endif
    enddo

   if(node_list%node(inode)%constrained==.true.) then
	     do j = 1, 2
		 index1 = node_list%node(inode)%parents(j)
		  do k=1, n_order+1

                       if ( (node_list%node(index1)%index(k) .ge. index_min(my_id+1)) .and. &
                        (node_list%node(index1)%index(k) .le. index_max(my_id+1)) ) then
                        ELM_is_local = .true.
                        goto 10
                      endif
                  enddo
	    end do     
	    	    
   end if

  enddo
  if (ELM_is_local) then
10   inext = inext + 1
     local_elms(inext) = i
  endif

enddo

n_local_ELMs = inext

!write(*,'(i4,A,20i8)') my_id,' n_local_elms  : ',n_local_elms,element_list%n_elements

return
end

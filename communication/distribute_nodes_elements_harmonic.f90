subroutine distribute_nodes_elements_harmonic(my_id,m_cpu,n_cpu,node_list,element_list, &
                                    local_elms, n_local_elms, index_min_harm, index_max_harm)
!---------------------------------------------------------------------------------------------
! subroutine divides the nodes (not their individual dof) over m_cpu equal parts
!            builds local_elms, contain all elements with at least one node with 
!            one index between index_min and index_max
!---------------------------------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: flux_list

integer, intent(out) :: local_elms(*), n_local_elms
integer, intent(in) :: my_id, m_cpu, n_cpu
integer :: index_total, inode
!integer, intent(in) :: index_min_harm(*), index_max_harm(*)
integer, intent(out) :: index_min_harm(*), index_max_harm(*)
integer :: index_part, inext, i,j, k, iv,index1

logical :: elm_is_local
!integer, dimension(node_list%n_nodes) :: active_node
!integer                               :: n_active_nodes
if (my_id .eq.0) then
  write(*,*) '************************************'
  write(*,*) '*  distributing nodes harmonic     *'
  write(*,*) '************************************'
endif



index_total = -1
do inode=1,node_list%n_nodes
  index_total = max(index_total,maxval(node_list%node(inode)%index))
enddo

#ifdef PSV
   index_min_harm(my_id+1) = 1
   index_max_harm(my_id+1) = index_total!ndof_glob/(n_tor*n_var)
#endif 


!stop
!write(*,*) ' n_elements  : ',my_id,element_list%n_elements
!write(*,*) ' n_nodes     : ',my_id,node_list%n_nodes
!write(*,*) ' index_total : ',my_id,index_total

!#ifdef PSV
index_min_harm(1:n_cpu) = 0
index_max_harm(1:n_cpu) = 0

!!----------------------------- must really take into account the number of
!!elements contributing to each node
if(MOD(my_id,m_cpu).eq.0) index_min_harm(my_id+1) = 1
do i=1,n_cpu
  index_max_harm(i) = ((MOD(i-1,m_cpu)+1) * index_total) / m_cpu
enddo
do i=2,n_cpu
   if(MOD(i-1,m_cpu).ne.0) then
    index_min_harm(i) = index_max_harm(i-1) + 1
   endif
enddo 
if(mod(my_id+1,m_cpu).eq.0) index_max_harm(my_id+1) = index_total
write(*,'(A,3i6)') ' index_min_harm,index_max_harm :',my_id,index_min_harm(my_id+1),index_max_harm(my_id+1)
!#endif 

!----------------------------- must really take into account the number of
!elements contributing to each node
#ifdef PSV
index_min_harm(1:n_cpu) = 0
index_max_harm(1:n_cpu) = 0
index_min_harm(1) = 1
do i=1,n_cpu
  index_max_harm(i) = (i * index_total) / n_cpu
enddo
do i=2,n_cpu
  index_min_harm(i) = index_max_harm(i-1) + 1
enddo
if (my_id .eq. n_cpu-1) index_max_harm(my_id+1) = index_total
#endif 



!if(mod(my_id+1,m_cpu).eq.0) index_max_harm(my_id+1) = index_total
!!write(*,'(A,3i6)') ' index_min_harm,index_max_harm :',my_id,index_min_harm(my_id+1),index_max_harm(my_id+1)
!if (my_id .eq. n_cpu-1) index_max_harm(my_id+1) = index_total
! 
!----------------------------------------------- find the elements that have a local node
inext = 0

do i = 1, element_list%n_elements

  ELM_is_local = .false.

  L_IV: do iv=1,n_vertex_max

    inode = element_list%element(i)%vertex(iv)

    do k=1, n_order+1

      if ( (node_list%node(inode)%index(k) .ge. index_min_harm(my_id+1)) .and. &
           (node_list%node(inode)%index(k) .le. index_max_harm(my_id+1)) ) then
        ELM_is_local = .true.
        exit L_IV
      endif
    enddo

    if(node_list%node(inode)%constrained) then
	     do j = 1, 2
		 index1 = node_list%node(inode)%parents(j)
		  do k=1, n_order+1

                       if ( (node_list%node(index1)%index(k) .ge. index_min_harm(my_id+1)) .and. &
                        (node_list%node(index1)%index(k) .le. index_max_harm(my_id+1)) ) then
                        ELM_is_local = .true.
                        exit L_IV
                      endif
                  enddo
	    end do     
	    	    
    end if

  end do L_IV
  if (ELM_is_local) then
     inext = inext + 1
     local_elms(inext) = i 
!     if (my_id.eq.0) write(*,*) local_elms(inext)
  endif

enddo

n_local_ELMs = inext

write(*,'(i4,A,20i8)') my_id,' n_local_elms  : ',n_local_elms,element_list%n_elements


return
end

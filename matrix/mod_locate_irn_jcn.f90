module mod_locate_irn_jcn
implicit none
contains
subroutine locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position)
!**************************************************************************
! subroutine finds the position in the global matrix of the index of      *
! node1 and node2 (this is the index per block)                           *
!                                                                         *
! search to be replaced by binary search                                  *
!**************************************************************************
use global_distributed_matrix
integer :: index_node1, index_node2, index_min, index_max, ijA_position, i, index1_local
logical :: found_index

found_index = .false.

index1_local = index_node1 - index_min + 1

!write(*,'(A,8i8)') ' LOCATE : ',index_node1,index_min,index_max,index1_local

do i=1,ijA_size(index1_local)           ! replace by binary search?

  if (irn_jcn(index1_local,i) .eq. index_node2) then
    ijA_position = ijA_index(index1_local,i)
    found_index = .true.
    exit
  endif

enddo

if (.not.found_index) then

  write(*,*) ' FATAL locate_irn_jcn : index not found ',index_node1,index_node2

  do i=1,ijA_size(index1_local)           ! replace by binary search?

    write(*,*) i,irn_jcn(index1_local,i)
 stop
  enddo

endif

return
end subroutine locate_irn_jcn
end module mod_locate_irn_jcn

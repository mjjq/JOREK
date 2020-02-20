module mod_locate_irn_jcn
implicit none
contains
subroutine locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position,& 
                          ijA_index, ijA_size, irn_jcn)
!**************************************************************************
! subroutine finds the position in the global matrix of the index of      *
! node1 and node2 (this is the index per block)                           *
!                                                                         *
! search to be replaced by binary search                                  *
!**************************************************************************
use mpi_mod
integer :: index_node1, index_node2, index_min, index_max, ijA_position, i, index1_local
logical :: found_index
integer :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
integer :: my_id, ierr


call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)


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

    write(*,*) 'my_id, i,irn_jcn(index1_local,i)', my_id, i,irn_jcn(index1_local,i)
 stop
  enddo

endif

return
end subroutine locate_irn_jcn
end module mod_locate_irn_jcn

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
!use global_distributed_matrix
use mpi_mod
integer :: index_node1, index_node2, index_min, index_max, ijA_position, i, index1_local
logical :: found_index
integer :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
integer :: my_id, rank, ierr


    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    my_id = rank


found_index = .false.

index1_local = index_node1 - index_min + 1

!write(*,'(A,8i8)') ' LOCATE : ',index_node1,index_min,index_max,index1_local

!!----add by PSV 
!if (ijA_size(index1_local) .eq. 0) ijA_size(index1_local) = 1
!!----add by PSV 

do i=1,ijA_size(index1_local)           ! replace by binary search?

  if (irn_jcn(index1_local,i) .eq. index_node2) then
    ijA_position = ijA_index(index1_local,i)
    found_index = .true.
    exit
  endif

enddo

if (.not.found_index) then

  !write(*,*) ' FATAL locate_irn_jcn : index not found ',index_node1,index_node2
  write(*,*) ' FATAL locate_irn_jcn : index not found '!,index_node1,index_node2
  write(*,*) 'my_id, index_node1, index_node2, index_min', my_id, index_node1,index_node2, index_min
  write(*,*) 'my_id, ijA_size(index1_local), irn_jcn(index1_local,i)', my_id, ijA_size(index1_local), irn_jcn(index1_local,1)  

  do i=1,ijA_size(index1_local)           ! replace by binary search?

    write(*,*) 'my_id, i,irn_jcn(index1_local,i)', my_id, i,irn_jcn(index1_local,i)
 stop
  enddo

endif

return
end subroutine locate_irn_jcn
end module mod_locate_irn_jcn

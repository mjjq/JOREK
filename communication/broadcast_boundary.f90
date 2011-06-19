subroutine Broadcast_boundary(my_id,boundary_list,bnd_node_list)
!----------------------------------------------------------
! subroutine to broadcast all the nodes and elements of the boundary
!----------------------------------------------------------
use tr_module
use data_structure

implicit none

include 'mpif.h'                                       ! MPI fortran include file

type (type_bnd_element_list) :: boundary_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element)      :: aboundary
type (type_bnd_node)         :: abnd_node
integer                      :: my_id, ife, ind, ierr, position, bufsize, IDBL_EXT, INT_EXT, ILOG_EXT
integer, allocatable         :: buffer(:)


call MPI_TYPE_EXTENT(MPI_INTEGER,INT_EXT,ierr)
call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,IDBL_EXT,ierr)

call MPI_BCAST(boundary_list%n_bnd_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(bnd_node_list%n_bnd_nodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

bufsize = boundary_list%n_bnd_elements * (10*INT_EXT + 4*IDBL_EXT) + &
          bnd_node_list%n_bnd_nodes    * (4*INT_EXT)

call tr_allocate(buffer,1,int(bufsize/INT_EXT),"bcastb_buffer")

if (my_id .eq. 0) then
  position = 0

  do ife=1,boundary_list%n_bnd_elements

    aboundary = boundary_list%bnd_element(ife)

    call MPI_PACK(aboundary%vertex,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%bnd_vertex,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%direction,4,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%element,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%side,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%size,4,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  enddo

  do ind=1,bnd_node_list%n_bnd_nodes

    abnd_node = bnd_node_list%bnd_node(ind)

    call MPI_PACK(abnd_node%index_jorek,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(abnd_node%index_starwall,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(abnd_node%direction,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  enddo
  
endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

  position = 0
  do ife=1,boundary_list%n_bnd_elements

    call MPI_UNPACK(buffer,bufsize,position,aboundary%vertex,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%bnd_vertex,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%direction,4,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%element,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%side,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%size,4,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    boundary_list%bnd_element(ife) = aboundary

  enddo

  do ind=1,bnd_node_list%n_bnd_nodes

    call MPI_UNPACK(buffer,bufsize,position,abnd_node%index_jorek,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,abnd_node%index_starwall,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,abnd_node%direction,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    
    bnd_node_list%bnd_node(ind) = abnd_node
  enddo
  
endif

call tr_deallocate(buffer,"bcastb_buffer")

return
end

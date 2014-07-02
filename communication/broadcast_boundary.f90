!> Broadcast all the nodes and elements of the boundary
subroutine broadcast_boundary(my_id,boundary_list,bnd_node_list)

use tr_module
use data_structure
use mpi_mod

implicit none


! --- Routine parameters
integer,                      intent(in)    :: my_id
type (type_bnd_element_list), intent(inout) :: boundary_list
type (type_bnd_node_list),    intent(inout) :: bnd_node_list

! --- Local variables
type (type_bnd_element) :: aboundary
type (type_bnd_node)    :: abnd_node
integer                 :: ife, ind, ierr, position, bufsize, IDBL_EXT, INT_EXT, ILOG_EXT
character, allocatable  :: buffer(:)

call MPI_BCAST(boundary_list%n_bnd_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(bnd_node_list%n_bnd_nodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! --- Prepare broadcast buffer
call MPI_PACK_SIZE(1,MPI_INTEGER,MPI_COMM_WORLD,INT_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IDBL_EXT,ierr)
bufsize = boundary_list%n_bnd_elements * (10*INT_EXT + 4*IDBL_EXT) + &
          bnd_node_list%n_bnd_nodes    * (6*INT_EXT)
allocate(buffer(bufsize))
call tr_register_mem(bufsize,"bcastb_buffer")

if (my_id == 0) then
  position = 0

  do ife=1,boundary_list%n_bnd_elements

    aboundary = boundary_list%bnd_element(ife)

    call MPI_PACK(aboundary%vertex,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%bnd_vertex,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%direction,4,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%element,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%side,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%size,4,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  end do

  do ind=1,bnd_node_list%n_bnd_nodes

    abnd_node = bnd_node_list%bnd_node(ind)

    call MPI_PACK(abnd_node%index_jorek,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(abnd_node%index_starwall,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(abnd_node%direction,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(abnd_node%n_dof,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  end do
  
end if

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id /= 0) then

  position = 0
  do ife=1,boundary_list%n_bnd_elements

    call MPI_UNPACK(buffer,bufsize,position,aboundary%vertex,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%bnd_vertex,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%direction,4,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%element,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%side,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%size,4,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    boundary_list%bnd_element(ife) = aboundary

  end do

  do ind=1,bnd_node_list%n_bnd_nodes

    call MPI_UNPACK(buffer,bufsize,position,abnd_node%index_jorek,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,abnd_node%index_starwall,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,abnd_node%direction,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,abnd_node%n_dof,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    
    bnd_node_list%bnd_node(ind) = abnd_node
  end do
  
end if

call tr_unregister_mem(bufsize,"bcastb_buffer")
deallocate(buffer)

return
end subroutine broadcast_boundary

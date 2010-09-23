subroutine Broadcast_boundary(my_id,boundary_list)
!----------------------------------------------------------
! subroutine to broadcast all the nodes in the point_list
!----------------------------------------------------------
use data_structure
implicit none

type (type_bnd_element_list) :: boundary_list
include 'mpif.h'                                       ! MPI fortran include file
type (type_bnd_element)      :: aboundary
integer                      :: my_id, ife, ierr, position, bufsize, IDBL_EXT, INT_EXT, ILOG_EXT
integer, allocatable         :: buffer(:)

!  type type_boundary                                  ! type definition for one boundary element (1D element)
!    integer :: vertex(2)                              ! the nodes of the corners
!    integer :: direction(2)                           ! indicates which direction of the nodes is along the boundary (2 or 3)
!    integer :: element                                ! boundary element is part of this element
!    integer :: side                                   ! boundary element corresponds to this side of the originating element
!    real*8  :: size(2,2)                              ! the size of the vectors at each vertex of the element
!  endtype type_boundary

!  type type_boundary_list                             ! type definition for a list of boundary elements
!    integer :: n_boundary                             ! the number of boundary elements in the list
!    type (type_boundary)  :: boundary(n_elements_max)  ! the list of boundary elements
!  endtype type_boundary_list


call MPI_TYPE_EXTENT(MPI_INTEGER,INT_EXT,ierr)
call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,IDBL_EXT,ierr)

call MPI_BCAST(boundary_list%n_bnd_elements,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

bufsize = boundary_list%n_bnd_elements * (8*INT_EXT + 4*IDBL_EXT)

if (allocated(buffer)) deallocate(buffer)
allocate(buffer(bufsize/ INT_EXT))

if (my_id .eq. 0) then
  position = 0

  do ife=1,boundary_list%n_bnd_elements

    aboundary = boundary_list%bnd_element(ife)

    call MPI_PACK(aboundary%vertex,2,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%direction,4,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%element,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%side,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(aboundary%size,4,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  enddo

endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

  position = 0
  do ife=1,boundary_list%n_bnd_elements

    call MPI_UNPACK(buffer,bufsize,position,aboundary%vertex,2,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%direction,4,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%element,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%side,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,aboundary%size,4,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    boundary_list%bnd_element(ife) = aboundary

  enddo

endif

deallocate(buffer)

return
end

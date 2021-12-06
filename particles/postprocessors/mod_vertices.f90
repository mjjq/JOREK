!> the module mod_vertices contains all datatypes
!> and procedures common to light and gather points
module mod_vertices
implicit none

private
public :: vertices

!> Variable and type definitions ------------------------------------------
!> vertices: abstract class containing the basic types
!> and procedures defining light or gather points
type,abstract :: vertices
  !> type variables
  integer,parameter     :: n_x=3             !< number of coordinates
  integer               :: n_vertices        !< total number of vertices
  integer               :: n_active_vertices !< number of active verticies
  integer               :: n_property_vertex !< total number of properties per vertex
  real*8,dimension(:,:) :: x                 !< position of the point (cartesian)
  real*8,dimension(:,:) :: properties        !< properties of the vertices
  contains
  !> type procedures
  procedure,pass(vertices) :: allocate_vertices
  procedure,pass(vertices) :: deallocate_vertices
  procedure,pass(vertices) :: resize_vertices_noloss
  procedure,pass(vertices) :: fit_tables_to_active_vertices
end type vertices
!> Interfaces -------------------------------------------------------------
contains

!> Procedures -------------------------------------------------------------
!> allocate vertex tables and initilise them to 0. Data loss is expected
!> for preallocated tables
!> inputs:
!>   vert_inout:        (vertices) vertices to be allocated
!>   n_vertices:        (integer) number of vertices
!> outputs:
!>   vert_inout: (vertices) allocated vertices
subroutine allocate_vertices(vert_inout,n_vertices)
  implicit none
  !> inputs
  integer,intent(in) :: n_vertices
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  !> check,allocate and initialize to 0 x and properties array
  if(allocated(vert_inout%x)) deallocate(vert_inout%x)
  if(allocated(vert_inout%properties)) deallocate(vert_inout%properties)
  allocate(vert_inout%x(vert_inout%n_x,n_vertices))
  allocate(vert_inout%properties(vert_inout%n_property_vertex,n_vertices))
  vert_inout%x = 0.d0; vert_inout%properties = 0.d0;
  vert_inout%n_vertices = n_vertices; vert_inout%n_active_vertices = 0;
end subroutine allocate_vertices

!> deallocate vertex tables. Data loss is expected
!> inputs:
!>   vert_inout: (vertices) allocated vertices
!> outputs:
!>   vert_inout: (vertices) deallocated vertices
subroutine deallocate_vertices(vert_inout)
  implicit none
  !> input-outputs
  class(vertices),intent(inout) :: vert_inout
  if(allocated(vert_inout%x)) deallocate(vert_inout%x)
  if(allocated(vert_inout%properties)) deallocate(vert_inout%properties)
  vert_inout%n_vertices = 0; vert_inout%n_active_vertices = 0;
end subroutine deallocate_vertices

!> resize the vertex tables without loss of data. This operation can be
!> very slow because it requires a copy of all tables data
!> inputs:
!>   vert_inout:   (vertices) vertex table with old sizes
!>   n_vertex_new: (integer) new size of vertex tables
!>   ifail:        (integer) failed to resize tables if =11
!> outputs:
!>  vert_inout:    (vertices) resized vertex tables
!>  ifail:         (integer) failed to resize if =11
subroutine resize_vertices_noloss(vert_inout,n_vertex_new,ifail)
  implicit none
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  integer,intent(inout)         :: ifail
  !> inputs
  integer,intent(in)            :: n_vertex_new 
  !> variables
  integer :: ii
  real*8,dimension(vert_inout%n_x,vert_inout%n_active_vertices) :: x_table
  real*8,dimension(vert_inout%n_property_vertex,vert_inout%n_active_vertices) :: property_table

  !> check for possible data losses
  if(n_vertex.lt.vert_inout%n_active_vertices) then
    write(*,*) "Try to resize vertex tables but possible data loss detected!"
    ifail = 11
    return
  endif
  !> copy data and resize tables
  !$omp parallel default(private) shared(vert_inou,n_vertex_new,&
  !$omp x_table,property_table)
  !$omp do
  do ii=1,vert_inout%n_active_vertices
    x_table(:,ii) = vert_inout%x(:,ii)
    property_table(:,ii) = vert_inout%properties(:,ii)
  enddo
  !$omp end do
  !$omp single
  call vert_inout%allocate_vertices(n_vertex_new)
  !$omp end single
  !$omp do
  do ii=1,vert_inout%n_active_vertices
    vert_inout%x(:,ii) = x_table(:,ii)
    vert_inout%properties(:,ii) = property_table(:,ii)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine resize_vertices_noloss

!> fit the vertices table to the number of active tables.
!> This operation can be very slow because it requires 
!> a copy of all tables data
!> inputs:
!>   vert_inout: (vertices) vertex table with empty vertices
!>   ifail:      (integer) failed to resize tables if =11
!> outputs:
!>  vert_inout:  (vertices) fitted vertex tables
!>  ifail:       (integer) failed to resize if =11
subroutine fit_tables_to_active_vertices(vert_inout,ifail)
  implicit none
  class(vertices),intent(inout) :: vert_inout
  integer,intent(inout)         :: ifail
  call vert_inout%resize_vertices_noloss(vert_inout%n_active_vertices,ifail)
end subroutine fit_tables_to_active_vertices

!>-------------------------------------------------------------------------
end module mod_vertices

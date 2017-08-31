!> Module for creating an RTree of elements and searching in this tree
module mod_element_rtree
use iso_c_binding
use data_structure
implicit none
public populate_element_rtree, nearby_elements, elements_containing_point, initialized

logical :: initialized = .false.

interface
  !> Name is element_rtree to match filename `element_rtree.cpp`.
  !> `void PopulateTree(int nelm, double minx[], double miny[], double maxx[], double maxy[])`
  subroutine element_rtree(n, minx, miny, maxx, maxy) bind(C,name="PopulateTree")
    import C_DOUBLE, C_INT
    integer(C_INT), intent(in), value :: n
    real(C_DOUBLE), intent(in), dimension(n)  :: minx, miny, maxx, maxy
  end subroutine element_rtree

  !> `int NumElementsInRect(double minx, double miny, double maxx, double maxy)`
  function num_elements_in_rect(minx, miny, maxx, maxy) bind(C,name='NumElementsInRect')
    import C_DOUBLE, C_INT
    real(C_DOUBLE), intent(in), value         :: minx, miny, maxx, maxy
    integer(C_INT) :: num_elements_in_rect
  end function num_elements_in_rect
  !> `int ElementsInRect(double minx, double miny, double maxx, double maxy, int *ielm)`
  function elements_in_rect(minx, miny, maxx, maxy, ielm) bind(C,name='ElementsInRect')
    import C_DOUBLE, C_INT
    real(C_DOUBLE), intent(in), value         :: minx, miny, maxx, maxy
    integer(C_INT), intent(out), dimension(*) :: ielm
    integer(C_INT) :: elements_in_rect
  end function elements_in_rect
end interface

contains

!> Populate the RTree with the squares containing elements
!> x=R, y=Z
!>
!> Only one call to this routine can be made simultaneously! It uses a static
!> data structure under the hood.
subroutine populate_element_rtree(node_list, element_list)
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  real(C_DOUBLE), dimension(:), allocatable :: minx, miny, maxx, maxy
  real*8 :: rmin, rmax, zmin, zmax
  integer :: i, n
  n = element_list%n_elements
  allocate(minx(n), miny(n), maxx(n), maxy(n))
  do i=1,n
    call RZ_minmax(node_list, element_list, i, rmin, rmax, zmin, zmax)
    minx(i) = real(rmin, kind=C_DOUBLE)
    maxx(i) = real(rmax, kind=C_DOUBLE)
    miny(i) = real(zmin, kind=C_DOUBLE)
    maxy(i) = real(zmax, kind=C_DOUBLE)
  end do
  write(*,*) "Initializing RTree"
  ! this cleans out the tree before insertion
  call element_rtree(int(n,C_INT), minx, miny, maxx, maxy)
  initialized = .true.
end subroutine populate_element_rtree

!> Find probable neighbours of element i and return their indices.
!> This is done by taking the bounding box of an element and expanding
!> it slightly (10^-6 in absolute value) and returning all elements in
!> this box, except element i
subroutine nearby_elements(node_list, element_list, i_elm, i_nearby)
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  integer, intent(in)                 :: i_elm
  integer, dimension(:), allocatable, intent(out) :: i_nearby
  integer, parameter :: n_nearby_C = 500
  integer(C_int), dimension(n_nearby_C) :: i_nearby_C ! large placeholder. Should be larger than n_tht (or n_pol)

  integer :: i
  integer, dimension(n_vertex_max) :: vertices
  real(C_DOUBLE) :: minx, miny, maxx, maxy
  integer(C_INT) :: num_elements

  if (.not. initialized) then
    write(*,*) "ERROR: trying to find nearby elements in uninitialized tree, returning nothing"
    num_elements = 0
  else
    vertices = element_list%element(i_elm)%vertex(:)
    minx = real(minval(node_list%node(vertices)%x(1,1)) - 1d-6, kind=C_DOUBLE)
    maxx = real(maxval(node_list%node(vertices)%x(1,1)) + 1d-6, kind=C_DOUBLE)
    miny = real(minval(node_list%node(vertices)%x(1,2)) - 1d-6, kind=C_DOUBLE)
    maxy = real(maxval(node_list%node(vertices)%x(1,2)) + 1d-6, kind=C_DOUBLE)

    num_elements = int(elements_in_rect(minx, miny, maxx, maxy, i_nearby_C))
  end if
  if (num_elements .gt. n_nearby_C) then
    write(*,*) "ERROR: Too many elements returned by elements_in_rect, expect memory corruption (rect=", minx,maxx, ',', miny,maxy, ") n=", num_elements
    return
  end if
  allocate(i_nearby(num_elements))
  i_nearby = int(i_nearby_C(1:num_elements))
end subroutine nearby_elements

!> Find elements that could probably contain this point.
subroutine elements_containing_point(R, Z, i_elms)
  real*8, intent(in)                              :: R, Z
  integer, dimension(:), allocatable, intent(out) :: i_elms

  real(C_DOUBLE) :: minx, miny, maxx, maxy
  integer(C_INT) :: num_elements
  integer, parameter :: n_nearby_C = 500
  integer(C_int), dimension(n_nearby_C) :: i_nearby_C ! large placeholder. Should be larger than n_tht (or n_pol)

  if (.not. initialized) then
    write(*,*) "ERROR: trying to find point in uninitialized tree, returning nothing"
    num_elements = 0
  else
    minx = real(R, kind=C_DOUBLE)
    maxx = real(R, kind=C_DOUBLE)
    miny = real(Z, kind=C_DOUBLE)
    maxy = real(Z, kind=C_DOUBLE)

    ! This calls the search routine twice! To get around that either pass a large enough array alway
    ! or implement it as a mask/bitfield of the total number of elements.
    ! I expect this to be unnecessary, as the speed improvement here is dramatic enough.
    num_elements = int(elements_in_rect(minx, miny, maxx, maxy, i_nearby_C))
  end if
  if (num_elements .gt. n_nearby_C) then
    write(*,*) "ERROR: Too many elements returned by elements_containing_point, expect memory corruption, point=(", R, ',', Z, ") n=", num_elements
    return
  end if
  allocate(i_elms(num_elements))
  i_elms = int(i_nearby_C(1:num_elements))
end subroutine elements_containing_point
end module mod_element_rtree

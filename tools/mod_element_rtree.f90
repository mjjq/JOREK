!> Module for creating an RTree of elements and searching in this tree
module mod_element_rtree
use iso_c_binding
use data_structure
implicit none
public

interface
  !> Name is element_rtree to match filename `element_rtree.cpp`.
  !> `void PopulateTree(int nelm, double minx[], double miny[], double maxx[], double maxy[])`
  subroutine element_rtree(n, minx, miny, maxx, maxy) bind(C,name="PopulateTree")
    import C_DOUBLE, C_INT
    integer(C_INT), intent(in), value :: n
    real(C_DOUBLE), intent(in), dimension(n)    :: minx, miny, maxx, maxy
  end subroutine element_rtree

  !> `int NumElementsInRect(double minx, double miny, double maxx, double maxy)`
  function num_elements_in_rect(minx, miny, maxx, maxy) bind(C,name='NumElementsInRect')
    import C_DOUBLE, C_INT
    real(C_DOUBLE), intent(in), value           :: minx, miny, maxx, maxy
    integer(C_INT) :: num_elements_in_rect
  end function num_elements_in_rect
  !> `void ElementsInRect(double minx, double miny, double maxx, double maxy, int *ielm)`
  subroutine elements_in_rect(minx, miny, maxx, maxy, ielm) bind(C,name='ElementsInRect')
    import C_DOUBLE, C_INT
    real(C_DOUBLE), intent(in), value           :: minx, miny, maxx, maxy
    integer(C_INT), intent(inout), dimension(*) :: ielm
  end subroutine elements_in_rect
end interface

contains

!> Populate the RTree with the squares containing nodes of elements.
!> x=R, y=Z
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
  call element_rtree(int(n,C_INT), minx, miny, maxx, maxy)
end subroutine populate_element_rtree

!> Find probable neighbours of element i and return their indices.
!> This is done by taking the bounding box of an element and expanding
!> it slightly (10^-6 in absolute value) and returning all elements in
!> this box, except element i
subroutine nearby_elements(node_list, element_list, i_elm, i_nearby)
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  integer, intent(in)                 :: i_elm
  integer(C_INT), dimension(:), allocatable, intent(out) :: i_nearby

  integer :: i
  integer, dimension(n_vertex_max) :: vertices
  real(C_DOUBLE) :: minx, miny, maxx, maxy
  integer(C_INT) :: num_elements

  vertices = element_list%element(i_elm)%vertex(:)
  minx = real(minval(node_list%node(vertices)%x(1,1)) - 1d-6, kind=C_DOUBLE)
  maxx = real(maxval(node_list%node(vertices)%x(1,1)) + 1d-6, kind=C_DOUBLE)
  miny = real(minval(node_list%node(vertices)%x(1,2)) - 1d-6, kind=C_DOUBLE)
  maxy = real(maxval(node_list%node(vertices)%x(1,2)) + 1d-6, kind=C_DOUBLE)

  ! This calls the search routine twice! To get around that either pass a large enough array alway
  ! or implement it as a mask/bitfield of the total number of elements.
  ! I expect this to be unnecessary, as the speed improvement here is dramatic enough.
  num_elements = num_elements_in_rect(minx, miny, maxx, maxy)
  allocate(i_nearby(num_elements))
  call elements_in_rect(minx, miny, maxx, maxy, i_nearby)
end subroutine nearby_elements

!> Find elements that could probably contain this point.
!> Do this by maxing a small box of 10^-6 around this point and finding all elements of which the RZ boundingbox overlaps.
!> This could be a few, but not more than n_pol (if exactly on axis)
subroutine elements_containing_point(R, Z, i_elms)
  real*8, intent(in)                  :: R, Z
  integer(C_INT), dimension(:), allocatable, intent(out) :: i_elms

  real(C_DOUBLE) :: minx, miny, maxx, maxy
  integer(C_INT) :: num_elements

  minx = real(R - 1d-6, kind=C_DOUBLE)
  maxx = real(R + 1d-6, kind=C_double)
  miny = real(Z - 1d-6, kind=C_double)
  maxy = real(Z + 1d-6, kind=C_double)

  ! This calls the search routine twice! To get around that either pass a large enough array alway
  ! or implement it as a mask/bitfield of the total number of elements.
  ! I expect this to be unnecessary, as the speed improvement here is dramatic enough.
  num_elements = num_elements_in_rect(minx, miny, maxx, maxy)
  allocate(i_elms(num_elements))
  call elements_in_rect(minx, miny, maxx, maxy, i_elms)
end subroutine elements_containing_point
end module mod_element_rtree

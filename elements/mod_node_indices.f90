!> With bi-cubic elements, ordering is easy, 1-4 nodes.
!> With n_order=5 and higher, ordering becomes quite complex
!> This module contains a set of routines to find the indices for you
module mod_node_indices

use mod_parameters

contains

!> calculate the node indices as a function of (s,t)-coord indices
subroutine calculate_node_indices(node_indices)

  implicit none

  integer, intent(inout) :: node_indices((n_order+1)/2,(n_order+1)/2)
  integer                :: k, l, i, j, counter

  node_indices = 0
  counter = 1
  ! --- We do it square by square
  do k = 1,(n_order+1)/2
    ! --- For each square, we do i row and j column
    do l = 1,k
        ! --- First the row
        i = k
        j = l
        node_indices(i,j) = counter
        counter = counter + 1
        ! --- Then the column
        i = l
        j = k
        if (i .eq. j) cycle ! don't record corners twice
        node_indices(i,j) = counter
        counter = counter + 1
    enddo
  enddo

  return
end subroutine calculate_node_indices


!> given the node index, find the (s,t)-coord indices of a node
subroutine get_node_coords_from_index(node_indices, index, i, j)

  implicit none

  integer, intent(in)    :: node_indices((n_order+1)/2,(n_order+1)/2), index
  integer, intent(inout) :: i, j
  integer                :: k, l
  logical                :: found

  i = 0 ; j = 0
  found = .false.
  do k = 1,(n_order+1)/2
    do l = 1,(n_order+1)/2
      if (node_indices(k,l) .eq. index) then
        i = k
        j = l
        found = .true.
        exit
      endif
    enddo
    if (found) exit
  enddo

  return
end subroutine get_node_coords_from_index





end module mod_node_indices

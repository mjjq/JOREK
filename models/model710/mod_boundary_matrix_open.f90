module mod_boundary_matrix_open
  implicit none
contains
subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(2)        ! the two nodes containing the boundary nodes

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS

integer    :: vertex(2), direction(2), xcase2
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
logical    :: xpoint2

return
end subroutine boundary_matrix_open
end module mod_boundary_matrix_open

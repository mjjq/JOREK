subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)
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

real*8     :: ELM(n_vertex_max*n_var*(n_order+1)*n_tor,n_vertex_max*n_var*(n_order+1)*n_tor)
real*8     :: RHS(n_vertex_max*n_var*(n_order+1)*n_tor)

integer    :: vertex(2), direction(2)
real*8     :: psi_axis, psi_bnd, Z_xpoint
logical    :: xpoint2

return
end

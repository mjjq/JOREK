module elm_fft
use parameters

  real*8 :: ELM_p(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: ELM_n(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: ELM_k(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: ELM_kn(n_plane,n_vertex_max*n_var*(n_order+1),n_vertex_max*n_var*(n_order+1)*n_tor)
  real*8 :: RHS_p(n_plane,n_vertex_max*n_var*(n_order+1))
  real*8 :: RHS_k(n_plane,n_vertex_max*n_var*(n_order+1))

!$OMP THREADPRIVATE(ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k)
end module

module mod_elt_matrix_fft
contains
subroutine element_matrix_fft(element, nodes, xpoint2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use elm_fft
implicit none
 
type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS
integer, intent(in) :: tid
real*8     :: psi_axis, psi_bnd, Z_xpoint
logical    :: xpoint2

return
end subroutine element_matrix_fft
end module mod_elt_matrix_fft
      

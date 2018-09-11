! This module contains nothing (just a wrapper) but it is needed by construct_matrix for the other models.
! Can be removed once the other models have also combined element_matrix and element_matrix_fft.
module mod_elt_matrix
contains

  subroutine element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
  !--------------------------------------------------------------------------
  ! This is just a wrapper to the real routine since I combined both into one
  !--------------------------------------------------------------------------

    use data_structure
    use mod_elt_matrix_fft

    implicit none

    type (type_element) 	      :: element
    type (type_node)		      :: nodes(n_vertex_max)

    integer    :: xcase2
    logical    :: xpoint2
    real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
    real*8, dimension (:,:), pointer  :: ELM
    real*8, dimension (:)  , pointer  :: RHS
    integer, intent(in) 	      :: tid

    call element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)

    return

  end subroutine element_matrix

end module mod_elt_matrix


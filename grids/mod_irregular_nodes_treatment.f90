!> Functionality for improving numerical stability at irregular nodes.
!! - Switched on by the logical namelist parameter treat_irregular_nodes.
!! - When used, the dofs of all irregular nodes are replaced by the following:
!!   1: value
!!   2: d/dR
!!   3: d/dZ
!!   4: d^2/dRdZ (to keep free)
!! - For further details, refer to the JOREK wiki at ??
module mod_irregular_nodes_treatment

contains

! Transform basis functions for the irregular nodes. This will solve for new degrees of freedom at the irregular nodes.
subroutine transform_basis_for_irregular_element_poisson(nodes, ELM, RHS, ivar_in, ivar_out, i_harm)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
real*8,          intent(inout) :: ELM(1:n_vertex_max*n_degrees, 1:n_vertex_max*n_degrees)
real*8,          intent(inout) :: RHS(1:n_vertex_max*n_degrees)
integer ,        intent(in)    :: ivar_in, ivar_out, i_harm
! --- routine parameters
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
real*8  :: Pmat_i(1:n_degrees, 1:n_degrees), Pmat_j(1:n_degrees, 1:n_degrees)
real*8  :: rhs_i(1:n_degrees), elm_ij(1:n_degrees, 1:n_degrees)

do iv = 1, n_vertex_max

  call construct_transformation_matrix(nodes(iv), Pmat_i)

  ! Extract RHS associated with a vertex iv: rhs_iv
  do io = 1, n_degrees
    index_iv_io = (iv-1)*n_degrees + io
    rhs_i(io)   = RHS(index_iv_io)
  enddo
  ! transform test functions for rhs_v
  rhs_i = matmul(transpose(Pmat_i), rhs_i)
  ! fill the updated entries in RHS vector
  do io = 1, n_degrees
    index_iv_io = (iv-1)*n_degrees + io    
    RHS(index_iv_io) = rhs_i(io)
  enddo

  do jv = 1, n_vertex_max

    call construct_transformation_matrix(nodes(jv), Pmat_j)

    ! Extract ELM associated with a vertex iv and jv: elm_iv_jv
    do io = 1, n_degrees
       index_iv_io = (iv-1)*n_degrees + io
       do jo = 1, n_degrees
          index_jv_jo = (jv-1)*n_degrees + jo
          elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
       enddo
    enddo   
    ! transform test and basis functions for elm_iv_jv 
    elm_ij = matmul(matmul(transpose(Pmat_i), elm_ij), Pmat_j)
    ! fill the updated entries in ELM matrix
    do io = 1, n_degrees
       index_iv_io = (iv-1)*n_degrees + io
       do jo = 1, n_degrees
          index_jv_jo = (jv-1)*n_degrees + jo
          ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
       enddo
    enddo

  enddo ! loop over jv

enddo ! loop over iv

end subroutine transform_basis_for_irregular_element_poisson

! Transform basis functions for the irregular nodes. This will solve for new degrees of freedom at the irregular nodes.
subroutine transform_basis_for_irregular_element(nodes, ELM, RHS, i_v, n_v, i_n, n_harm)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
real*8   :: ELM(1:n_tor*n_vertex_max*n_degrees*n_var, 1:n_tor*n_vertex_max*n_degrees*n_var)
real*8   :: RHS(1:n_tor*n_vertex_max*n_degrees*n_var)
integer,         intent(in)    :: n_v, i_v(n_v), n_harm, i_n(n_harm)

! --- routine parameters
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
integer :: ivar, jvar, im, in, n_tor_local
real*8  :: Pmat_i(1:n_degrees, 1:n_degrees), Pmat_j(1:n_degrees, 1:n_degrees)
real*8  :: rhs_i(1:n_degrees), elm_ij(1:n_degrees, 1:n_degrees)

n_tor_local = i_n(n_harm) - i_n(1) +1
  
do iv = 1, n_vertex_max

  call construct_transformation_matrix(nodes(iv), Pmat_i)
  
  do ivar = 1, n_v
  do im   = 1, n_harm
    ! Extract RHS associated with a vertex iv: rhs_iv
    do io = 1, n_degrees
      index_iv_io = n_tor_local*n_var*n_degrees*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local
      rhs_i(io)   = RHS(index_iv_io)
    enddo
    ! transform test functions for rhs_v
    rhs_i = matmul(transpose(Pmat_i), rhs_i)
    ! fill the updated entries in RHS vector
    do io = 1, n_degrees
      index_iv_io = n_tor_local*n_var*n_degrees*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local
      RHS(index_iv_io) = rhs_i(io)
    enddo
  enddo ! loop over im
  enddo ! loop over ivar

  do jv = 1, n_vertex_max

    call construct_transformation_matrix(nodes(jv), Pmat_j)

    ! Extract ELM associated with a vertex iv and jv: elm_iv_jv
    do ivar = 1, n_var
    do im   = 1, n_harm
    
    do jvar = 1, n_var
    do in   = 1, n_harm
    
       do io = 1, n_degrees
          index_iv_io = n_tor_local*n_var*n_degrees*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local
          do jo = 1, n_degrees
             index_jv_jo = n_tor_local*n_var*n_degrees*(jv-1) + n_tor_local*n_var*(jo-1) + (i_n(in)-i_n(1)+1) + (i_v(jvar)-1)*n_tor_local
             elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
          enddo
       enddo
       ! transform test and basis functions for elm_iv_jv
       elm_ij = matmul(matmul(transpose(Pmat_i), elm_ij), Pmat_j)
       ! fill the updated entries in ELM matrix
       do io = 1, n_degrees
          index_iv_io = n_tor_local*n_var*n_degrees*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local       
          do jo = 1, n_degrees
             index_jv_jo = n_tor_local*n_var*n_degrees*(jv-1) + n_tor_local*n_var*(jo-1) + (i_n(in)-i_n(1)+1) + (i_v(jvar)-1)*n_tor_local
             ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
          enddo
       enddo

    enddo ! loop over in
    enddo ! loop over jvar
       
    enddo ! loop over im
    enddo ! loop over ivar

  enddo ! loop over jv

enddo ! loop over iv

end subroutine transform_basis_for_irregular_element

! Since the irregular node treatment solves for new degrees of freedom, we need to
! transforms degrees of freedom to old ones.
subroutine new_to_old_dofs_on_irregular_node(node_list, vg, new_dofs, old_dofs)
use data_structure
use phys_module
implicit none
type(type_node_list),    intent(inout) :: node_list
integer,  intent(in) :: vg
real*8,  intent(in) :: new_dofs(1:n_degrees)
real*8,  intent(out):: old_dofs(1:n_degrees)
real*8  :: Pmat(1:n_degrees, 1:n_degrees), tmp(1:n_degrees)

call construct_transformation_matrix(node_list%node(vg), Pmat)
old_dofs = matmul(Pmat, new_dofs)
end subroutine new_to_old_dofs_on_irregular_node

subroutine construct_transformation_matrix(node, Pmat)
use data_structure
use phys_module
implicit none
type(type_node), intent(inout) :: node
real*8,          intent(out)   :: Pmat(1:n_degrees, 1:n_degrees)
integer :: j

Pmat = 0.d0
do j = 1, n_degrees
  Pmat(j,j) = 1.d0
enddo

if(n_order == 3)then

  if(node%irregular_node == 1)then
    Pmat(2,2) = node%x(1,2,1)
    Pmat(2,3) = node%x(1,2,2)

    Pmat(3,2) = node%x(1,3,1)
    Pmat(3,3) = node%x(1,3,2)

    Pmat(4,2) = node%x(1,4,1)
    Pmat(4,3) = node%x(1,4,2)
  elseif(node%irregular_node == 2)then
    Pmat(2,2) = node%x(1,2,1)
    Pmat(2,3) = node%x(1,2,2)

    Pmat(3,2) = node%x(1,3,1)
    Pmat(3,3) = node%x(1,3,2)

    Pmat(4,2) = node%x(1,4,1)
    Pmat(4,3) = node%x(1,4,2)    
  endif

elseif(n_order == 5)then

  Pmat(2,2) = node%x(1,2,1)
  Pmat(2,3) = node%x(1,2,2)

  Pmat(3,2) = node%x(1,3,1)
  Pmat(3,3) = node%x(1,3,2)

  Pmat(4,2) = node%x(1,4,1)
  Pmat(4,3) = node%x(1,4,2)
  Pmat(4,4) = node%x(1,2,1) * node%x(1,3,1) 
  Pmat(4,5) = node%x(1,2,1) * node%x(1,3,2) + node%x(1,3,1) * node%x(1,2,2)
  Pmat(4,6) = node%x(1,2,2) * node%x(1,3,2)
 
  Pmat(5,2) = node%x(1,5,1) 
  Pmat(5,3) = node%x(1,5,2)
  Pmat(5,4) = node%x(1,2,1) * node%x(1,2,1)
  Pmat(5,5) = node%x(1,2,1) * node%x(1,2,2) + node%x(1,2,1) * node%x(1,2,2)
  Pmat(5,6) = node%x(1,2,2) * node%x(1,2,2)

  Pmat(6,2) = node%x(1,6,1)
  Pmat(6,3) = node%x(1,6,2)
  Pmat(6,4) = node%x(1,3,1) * node%x(1,3,1)
  Pmat(6,5) = node%x(1,3,1) * node%x(1,3,2) + node%x(1,3,1) * node%x(1,3,2)
  Pmat(6,6) = node%x(1,3,2) * node%x(1,3,2)

  Pmat(7,2) = node%x(1,7,1)
  Pmat(7,3) = node%x(1,7,2)
  Pmat(7,4) = node%x(1,4,2)
  Pmat(7,5) = node%x(1,4,2)
  Pmat(7,6) = node%x(1,4,2)

  Pmat(8,2) = node%x(1,8,1)
  Pmat(8,3) = node%x(1,8,2)
  Pmat(8,4) = node%x(1,4,2)
  Pmat(8,5) = node%x(1,4,2)
  Pmat(8,6) = node%x(1,4,2)

  Pmat(9,4) = node%x(1,5,1) * node%x(1,6,1)
  Pmat(9,5) = node%x(1,5,1) * node%x(1,6,2) + node%x(1,6,1) * node%x(1,5,2)
  Pmat(9,6) = node%x(1,5,2) * node%x(1,6,2)
  
else

  write(*,*) 'Irregular point treatment is implemented for n_order 3 or 5 only'
  stop

endif  
        
end subroutine construct_transformation_matrix

end module mod_irregular_nodes_treatment

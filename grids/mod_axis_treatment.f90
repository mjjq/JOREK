module mod_axis_treatment

contains

! Transform basis functions for the axis nodes. This will solve for new degrees of freedom at the axis.
subroutine transform_basis_for_axis_element_poisson(nodes, ELM, RHS, ivar_in, ivar_out, i_harm)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
real*8,          intent(inout) :: ELM(1:n_vertex_max*(n_order+1), 1:n_vertex_max*(n_order+1))
real*8,          intent(inout) :: RHS(1:n_vertex_max*(n_order+1))
integer ,        intent(in)    :: ivar_in, ivar_out, i_harm
! --- routine parameters
integer :: axis_vertex1, axis_vertex4, dof1, dof2, dof3, dof4
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
real*8  :: Ptrans(1:4, 1:4), Pmat(1:4, 1:4), rhs_i(1:4), elm_ij(1:4, 1:4)

axis_vertex1 = 1 ; axis_vertex4 = 4
dof1 = 1 ; dof2 = 2 ; dof3 = 3 ; dof4 = 4

do iv = 1, n_vertex_max

  ! initialize transpose of the transformation matrix to identity
  Ptrans = 0.d0
  do io = 1, n_order+1
    Ptrans(io, io) = 1.d0
  enddo

  ! determine transpose of the transformation matrix
  if(iv==axis_vertex1 .or. iv==axis_vertex4)then
    Ptrans(dof2,dof2) = nodes(iv)%x(1,dof2,1) ; Ptrans(dof2,dof4) = nodes(iv)%x(1,dof4,1)
    Ptrans(dof3,dof3) = 0.d0
    Ptrans(dof4,dof2) = nodes(iv)%x(1,dof2,2) ; Ptrans(dof4,dof4) = nodes(iv)%x(1,dof4,2)
  endif

  ! Extract RHS associated with a vertex iv: rhs_iv
  do io = 1, n_order+1
    index_iv_io = (iv-1)*(n_order+1) + io
    rhs_i(io)   = RHS(index_iv_io)
  enddo
  ! transform test functions for rhs_v
  rhs_i = matmul( Ptrans, rhs_i)
  ! fill the updated entries in RHS vector
  do io = 1, n_order+1
    index_iv_io = (iv-1)*(n_order+1) + io    
    RHS(index_iv_io) = rhs_i(io)
  enddo

  do jv = 1, n_vertex_max

    ! initialize the transformation matrix to identity  
    Pmat = 0.d0
    do jo = 1, n_order+1
      Pmat(jo, jo) = 1.d0
    enddo

    ! determine the transformation matrix    
    if(jv==axis_vertex1 .or. jv==axis_vertex4)then
      Pmat(dof2,dof2)  = nodes(jv)%x(1,dof2,1) ; Pmat(dof2,dof4) = nodes(jv)%x(1,dof2,2)
      Pmat(dof3,dof3)  = 0.d0
      Pmat(dof4,dof2)  = nodes(jv)%x(1,dof4,1) ; Pmat(dof4,dof4) = nodes(jv)%x(1,dof4,2)
    endif

    ! Extract ELM associated with a vertex iv and jv: elm_iv_jv
    do io = 1, n_order+1
       index_iv_io = (iv-1)*(n_order+1) + io
       do jo = 1, n_order+1
          index_jv_jo = (jv-1)*(n_order+1) + jo
          elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
       enddo
    enddo   
    ! transform test and basis functions for elm_iv_jv 
    elm_ij = matmul(matmul(PTrans, elm_ij), Pmat)
    ! fill the updated entries in ELM matrix
    do io = 1, n_order+1
       index_iv_io = (iv-1)*(n_order+1) + io
       do jo = 1, n_order+1
          index_jv_jo = (jv-1)*(n_order+1) + jo
          ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
       enddo
    enddo

  enddo ! loop over jv

enddo ! loop over iv

end subroutine transform_basis_for_axis_element_poisson

! Transform basis functions for the axis nodes. This will solve for new degrees of freedom at the axis.
subroutine transform_basis_for_axis_element(nodes, ELM, RHS, i_v, n_v, i_n, n_harm)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
real*8   :: ELM(1:n_tor*n_vertex_max*(n_order+1)*n_var, 1:n_tor*n_vertex_max*(n_order+1)*n_var)
real*8   :: RHS(1:n_tor*n_vertex_max*(n_order+1)*n_var)
integer,         intent(in)    :: n_v, i_v(n_v), n_harm, i_n(n_harm)

! --- routine parameters
integer :: axis_vertex1, axis_vertex4, dof1, dof2, dof3, dof4
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
integer :: ivar, jvar, im, in, n_tor_local
real*8  :: Ptrans(1:4, 1:4), Pmat(1:4, 1:4), rhs_i(1:4), elm_ij(1:4, 1:4)

axis_vertex1 = 1 ; axis_vertex4 = 4
dof1 = 1 ; dof2 = 2 ; dof3 = 3 ; dof4 = 4

n_tor_local = i_n(n_harm) - i_n(1) +1

do iv = 1, n_vertex_max

  ! initialize transpose of the transformation matrix to identity
  Ptrans = 0.d0
  do io = 1, n_order+1
    Ptrans(io, io) = 1.d0
  enddo

  ! determine transpose of the transformation matrix
  if(iv==axis_vertex1 .or. iv==axis_vertex4)then
    Ptrans(dof2,dof2) = nodes(iv)%x(1,dof2,1) ; Ptrans(dof2,dof4) = nodes(iv)%x(1,dof4,1)
    Ptrans(dof3,dof3) = 0.d0    
    Ptrans(dof4,dof2) = nodes(iv)%x(1,dof2,2) ; Ptrans(dof4,dof4) = nodes(iv)%x(1,dof4,2)
  endif

  do ivar = 1, n_v
  do im   = 1, n_harm
    ! Extract RHS associated with a vertex iv: rhs_iv
    do io = 1, n_order+1
      index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local
      rhs_i(io)   = RHS(index_iv_io)
    enddo
    ! transform test functions for rhs_v
    rhs_i = matmul( Ptrans, rhs_i)
    ! fill the updated entries in RHS vector
    do io = 1, n_order+1
      index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local
      RHS(index_iv_io) = rhs_i(io)
    enddo
  enddo ! loop over im
  enddo ! loop over ivar

  do jv = 1, n_vertex_max

    ! initialize the transformation matrix to identity  
    Pmat = 0.d0
    do jo = 1, n_order+1
      Pmat(jo, jo) = 1.d0
    enddo

    ! determine the transformation matrix    
    if(jv==axis_vertex1 .or. jv==axis_vertex4)then
      Pmat(dof2,dof2)  = nodes(jv)%x(1,dof2,1) ; Pmat(dof2,dof4) = nodes(jv)%x(1,dof2,2)
      Pmat(dof3,dof3)  = 0.d0      
      Pmat(dof4,dof2)  = nodes(jv)%x(1,dof4,1) ; Pmat(dof4,dof4) = nodes(jv)%x(1,dof4,2)
    endif

    ! Extract ELM associated with a vertex iv and jv: elm_iv_jv
    do ivar = 1, n_var
    do im   = 1, n_harm
    
    do jvar = 1, n_var
    do in   = 1, n_harm
    
       do io = 1, n_order+1
          index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local
          do jo = 1, n_order+1
             index_jv_jo = n_tor_local*n_var*(n_order+1)*(jv-1) + n_tor_local*n_var*(jo-1) + (i_n(in)-i_n(1)+1) + (i_v(jvar)-1)*n_tor_local
             elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
          enddo
       enddo
       ! transform test and basis functions for elm_iv_jv
       elm_ij = matmul(matmul(PTrans, elm_ij), Pmat)
       ! fill the updated entries in ELM matrix
       do io = 1, n_order+1
          index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (i_n(im)-i_n(1)+1) + (i_v(ivar)-1)*n_tor_local       
          do jo = 1, n_order+1
             index_jv_jo = n_tor_local*n_var*(n_order+1)*(jv-1) + n_tor_local*n_var*(jo-1) + (i_n(in)-i_n(1)+1) + (i_v(jvar)-1)*n_tor_local
             ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
          enddo
       enddo

    enddo ! loop over in
    enddo ! loop over jvar
       
    enddo ! loop over im
    enddo ! loop over ivar

  enddo ! loop over jv

enddo ! loop over iv

end subroutine transform_basis_for_axis_element

! Since the axis treatment solves for new degrees of freedom, we need to
! transforms degrees of freedom to old ones.
subroutine new_to_old_dofs_poisson_on_the_axis(node_list, RHS)
use data_structure
use phys_module
implicit none
type(type_node_list),    intent(in) :: node_list
real*8,            intent(inout)    :: RHS(*)
integer :: i, k, in
integer :: dof2, dof4
integer :: index2, index4
real*8  :: Pmat(2,2), vec(2)

dof2 = 2 ; dof4 = 4

do i = 1, node_list%n_nodes
  if ( node_list%node(i)%axis_node ) then

    Pmat(1,1) = node_list%node(i)%x(1,dof2,1)  ;  Pmat(1,2) = node_list%node(i)%x(1,dof2,2)
    Pmat(2,1) = node_list%node(i)%x(1,dof4,1)  ;  Pmat(2,2) = node_list%node(i)%x(1,dof4,2)

    index2 = node_list%node(i)%index(dof2)
    index4 = node_list%node(i)%index(dof4)

    vec(1) = RHS(index2)
    vec(2) = RHS(index4)
    vec    = matmul(Pmat, vec)
    RHS(index2) = vec(1)
    RHS(index4) = vec(2)

  endif
enddo
end subroutine new_to_old_dofs_poisson_on_the_axis

! Since the axis treatment solves for new degrees of freedom, we need to
! transforms degrees of freedom to old ones.
subroutine new_to_old_dofs_on_the_axis(node_list, RHS)
use data_structure
use phys_module
implicit none
type(type_node_list),    intent(in) :: node_list
real*8,            intent(inout)    :: RHS(*)
integer :: i, k, in
integer :: dof2, dof4
integer :: index_node2, index_node4, index2, index4
real*8  :: Pmat(2,2), vec(2)

dof2 = 2 ; dof4 = 4
do i = 1, node_list%n_nodes
  if ( node_list%node(i)%axis_node ) then
    Pmat(1,1) = node_list%node(i)%x(1,dof2,1)  ;  Pmat(1,2) = node_list%node(i)%x(1,dof2,2)
    Pmat(2,1) = node_list%node(i)%x(1,dof4,1)  ;  Pmat(2,2) = node_list%node(i)%x(1,dof4,2)
    !do j=1,n_order+1
    index_node2 = node_list%node(i)%index(dof2)
    index_node4 = node_list%node(i)%index(dof4)    
    do k=1,n_var
      do in=1,n_tor
        index2 = n_tor*n_var * (index_node2 - 1) + n_tor*(k-1) + in
        index4 = n_tor*n_var * (index_node4 - 1) + n_tor*(k-1) + in         
        vec(1) = RHS(index2)
        vec(2) = RHS(index4)
        vec    = matmul(Pmat, vec)
        RHS(index2) = vec(1)
        RHS(index4) = vec(2)
      enddo
    enddo
    !enddo
  endif
enddo
end subroutine new_to_old_dofs_on_the_axis


end module mod_axis_treatment

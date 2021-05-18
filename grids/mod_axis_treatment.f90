module mod_axis_treatment

contains

! This subroutine transforms DoFs for a grid axis-node. The transfromation is from 
! new Dofs to old DoFs. This is done so that we can use old element matrix subroutines. 
pure subroutine transform_dofs_for_axis_node(node, i_v, n_v, i_n, n_harm, transform_deltas)
use data_structure
implicit none
type(type_node), intent(inout) :: node
integer,         intent(in)    :: n_v, i_v(n_v), n_harm, i_n(n_harm)
logical,         intent(in)    :: transform_deltas 
! --- routine parameters
integer :: iv, vg, ivar, in, dof2, dof4
real*8  :: Pmat(2,2), vec(2)

dof2 = 2 ; dof4 = 4

Pmat(1,1) = node%x(1,dof2,1)  ;  Pmat(1,2) = node%x(1,dof2,2)
Pmat(2,1) = node%x(1,dof4,1)  ;  Pmat(2,2) = node%x(1,dof4,2)

#if altcs
vec(1) = node%psi_eq(dof2)
vec(2) = node%psi_eq(dof4)
vec    = matmul(Pmat, vec)
node%psi_eq(dof2) = vec(1)
node%psi_eq(dof4) = vec(2)
#endif

do ivar = 1, n_v
   if(i_v(ivar) == 710)then
     vec(1) = node%Fprof_eq(dof2)
     vec(2) = node%Fprof_eq(dof4)
     vec    = matmul(Pmat, vec)
     node%Fprof_eq(dof2) = vec(1)
     node%Fprof_eq(dof4) = vec(2)
   elseif(i_v(ivar) == 711)then
     vec(1) = node%psi_eq(dof2)
     vec(2) = node%psi_eq(dof4)
     vec    = matmul(Pmat, vec)
     node%psi_eq(dof2) = vec(1)
     node%psi_eq(dof4) = vec(2)
   else
   do in   = 1, n_harm
      vec(1) = node%values(i_n(in),dof2,i_v(ivar))
      vec(2) = node%values(i_n(in),dof4,i_v(ivar))
      vec    = matmul(Pmat, vec)
      node%values(i_n(in),dof2,i_v(ivar)) = vec(1)
      node%values(i_n(in),dof2,i_v(ivar)) = vec(2)
    enddo       
   endif
enddo

if (transform_deltas) then
   do ivar = 1, n_v
   do in   = 1, n_harm
      vec(1) = node%deltas(i_n(in),dof2,i_v(ivar))
      vec(2) = node%deltas(i_n(in),dof4,i_v(ivar))
      vec    = matmul(Pmat, vec)
      node%deltas(i_n(in),dof2,i_v(ivar)) = vec(1)
      node%deltas(i_n(in),dof2,i_v(ivar)) = vec(2)
   enddo
   enddo
endif

end subroutine transform_dofs_for_axis_node

! This subroutine transforms basis functions in RHS and ELM for grid-axis-elements.
! This is done to achieve C1 continuity in (R, Z)-plane.
subroutine transform_basis_for_axis_element_poisson(nodes, ELM, RHS, ivar_in, ivar_out, i_harm)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
#define DIM0 n_vertex_max*(n_order+1)
real*8,          intent(inout) :: ELM(1:DIM0, 1:DIM0)
real*8,          intent(inout) :: RHS(1:DIM0)
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

! This subroutine transforms basis functions in RHS and ELM for grid-axis-elements.
! This is done to achieve C1 continuity in (R, Z)-plane.
subroutine transform_basis_for_axis_element(nodes, ELM, RHS, i_tor_min, i_tor_max)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
#define DIM0 n_tor*n_vertex_max*(n_order+1)*n_var
real*8,          intent(inout) :: ELM(1:DIM0, 1:DIM0)
real*8,          intent(inout) :: RHS(1:DIM0)
integer ,        intent(in)    ::  i_tor_min, i_tor_max
! --- routine parameters
integer :: axis_vertex1, axis_vertex4, dof1, dof2, dof3, dof4
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
integer :: ivar, jvar, im, in, n_tor_start, n_tor_end, n_tor_local
real*8  :: Ptrans(1:4, 1:4), Pmat(1:4, 1:4), rhs_i(1:4), elm_ij(1:4, 1:4)

axis_vertex1 = 1 ; axis_vertex4 = 4
dof1 = 1 ; dof2 = 2 ; dof3 = 3 ; dof4 = 4

n_tor_local = n_tor_end - n_tor_start +1

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

  do ivar = 1, n_var
  do im   = n_tor_start, n_tor_end
    ! Extract RHS associated with a vertex iv: rhs_iv
    do io = 1, n_order+1
      index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (im-n_tor_start+1) + (ivar-1)*n_tor_local
      rhs_i(io)   = RHS(index_iv_io)
    enddo
    ! transform test functions for rhs_v
    rhs_i = matmul( Ptrans, rhs_i)
    ! fill the updated entries in RHS vector
    do io = 1, n_order+1
      index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (im-n_tor_start+1) + (ivar-1)*n_tor_local
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
    do im   = n_tor_start, n_tor_end
    
    do jvar = 1, n_var
    do in   = n_tor_start, n_tor_end
    
       do io = 1, n_order+1
          index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (im-n_tor_start+1) + (ivar-1)*n_tor_local
          do jo = 1, n_order+1
             index_jv_jo = n_tor_local*n_var*(n_order+1)*(jv-1) + n_tor_local*n_var*(jo-1) + (in-n_tor_start+1) + (jvar-1)*n_tor_local
             elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
          enddo
       enddo
       ! transform test and basis functions for elm_iv_jv
       elm_ij = matmul(matmul(PTrans, elm_ij), Pmat)
       ! fill the updated entries in ELM matrix
       do io = 1, n_order+1
          index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + (im-n_tor_start+1) + (ivar-1)*n_tor_local
          do jo = 1, n_order+1
             index_jv_jo = n_tor_local*n_var*(n_order+1)*(jv-1) + n_tor_local*n_var*(jo-1) + (in-n_tor_start+1) + (jvar-1)*n_tor_local
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

! This subroutine is not used so far. We may not need it. Instead we can use transform_dofs_for_axis_node().
! This subroutine transforms DoFs for the nodes in the grid-axis-element.
!subroutine transform_dofs_for_axis_element(node_list, element, nodes, i_v, n_v, i_n, n_harm, transform_deltas)
!use data_structure
!implicit none
!type(type_node_list),  intent(in)    :: node_list
!type(type_element),    intent(inout) :: element
!type(type_node),       intent(inout) :: nodes(n_vertex_max)
!integer,               intent(in)    :: n_v, i_v(n_v), n_harm, i_n(n_harm)
!logical,               intent(in)    :: transform_deltas
!! --- routine parameters
!integer :: iv, vg, ivar, in, dof2, dof4
!real*8  :: Pmat(2,2), vec(2)
!
!dof2 = 2 ; dof4 = 4
!
!do iv = 1, 4, 3  ! loop over axis nodes only, 1 and 4 are axis nodes
!
!   vg = element%vertex(iv) ! global index of the vertex
!
!   Pmat(1,1) = nodes(iv)%x(1,dof2,1)  ;  Pmat(1,2) = nodes(iv)%x(1,dof2,2)
!   Pmat(2,1) = nodes(iv)%x(1,dof4,1)  ;  Pmat(2,2) = nodes(iv)%x(1,dof4,2)
!   
!#if fullmhd
!
!   vec(1) = node_list%node(vg)%Fprof_eq(dof2)
!   vec(2) = node_list%node(vg)%Fprof_eq(dof4)
!   vec    = matmul(Pmat, vec)
!   nodes(iv)%Fprof_eq(dof2) = vec(1)
!   nodes(iv)%Fprof_eq(dof4) = vec(2)
!
!   vec(1) = node_list%node(vg)%psi_eq(dof2)
!   vec(2) = node_list%node(vg)%psi_eq(dof4)
!   vec    = matmul(Pmat, vec)
!   nodes(iv)%psi_eq(dof2) = vec(1)
!   nodes(iv)%psi_eq(dof4) = vec(2)
!  
!#elif altcs
!
!   vec(1) = node_list%node(vg)%psi_eq(dof2)
!   vec(2) = node_list%node(vg)%psi_eq(dof4)
!   vec    = matmul(Pmat, vec)
!   nodes(iv)%psi_eq(dof2) = vec(1)
!   nodes(iv)%psi_eq(dof4) = vec(2)
!
!#endif
!
!   do ivar = 1, n_v
!   do in   = 1, n_harm
!      vec(1) = node_list%node(vg)%values(i_n(in),dof2,i_v(ivar))
!      vec(1) = node_list%node(vg)%values(i_n(in),dof4,i_v(ivar))
!      vec    = matmul(Pmat, vec)
!      nodes(iv)%values(in,dof2,ivar) = vec(1)
!      nodes(iv)%values(in,dof4,ivar) = vec(2)
!   enddo
!   enddo
!
!   if (transform_deltas) then
!   do ivar = 1, n_v
!   do in   = 1, n_harm
!      vec(1) = node_list%node(vg)%deltas(i_n(in),dof2,i_v(ivar))
!      vec(1) = node_list%node(vg)%deltas(i_n(in),dof4,i_v(ivar))
!      vec    = matmul(Pmat, vec)
!      nodes(iv)%deltas(in,dof2,ivar) = vec(1)
!      nodes(iv)%deltas(in,dof4,ivar) = vec(2)
!   enddo
!   enddo
!   endif
!
!enddo
!end subroutine transform_dofs_for_axis_element

! Following two subroutines 'transform_nodelist' and 'transform_back_nodelist' 
! transforms takes global node_list and transfrom axis nodes to and fro. Their
! use would make implementation of axis treament very easy. However, in the
! grid_xpoint subroutine 3rd and 4rth DoF for grid is set to zero: [x(1,3:4,:) = 0].
! Because of this back transformation is not invertible.

!subroutine transform_nodelist(node_list, i_tor_min, i_tor_max)
!use data_structure
!use phys_module
!implicit none
!integer :: vg, dof2, dof4, in, ivar
!type(type_node_list),  intent(inout) :: node_list
!integer ,              intent(in)    :: i_tor_min, i_tor_max
!real*8  :: Pmat(2,2), vec(2)
!
!dof2 = 2 ; dof4 = 4
!
!do vg = 1, node_list%n_nodes
!
!  if ( node_list%node(vg)%axis_node ) then
!
!    Pmat(1,1) = node_list%node(vg)%x(1,dof2,1)  ;  Pmat(1,2) = node_list%node(vg)%x(1,dof2,2)
!    Pmat(2,1) = node_list%node(vg)%x(1,dof4,1)  ;  Pmat(2,2) = node_list%node(vg)%x(1,dof4,2)
!
!!#if fullmhd
!!    vec(1) = node_list%node(vg)%Fprof_eq(dof2)
!!    vec(2) = node_list%node(vg)%Fprof_eq(dof4)
!!    vec    = matmul(Pmat, vec)
!!    node_list%node(vg)%Fprof_eq(dof2) = vec(1)
!!    node_list%node(vg)%Fprof_eq(dof4) = vec(2) 
!!
!!    vec(1) = node_list%node(vg)%psi_eq(dof2)
!!    vec(2) = node_list%node(vg)%psi_eq(dof4)
!!    vec    = matmul(Pmat, vec)
!!    node_list%node(vg)%psi_eq(dof2) = vec(1)
!!    node_list%node(vg)%psi_eq(dof4) = vec(2)
!!#elif altcs
!!    vec(1) = node_list%node(vg)%psi_eq(dof2)
!!    vec(2) = node_list%node(vg)%psi_eq(dof4)
!!    vec    = matmul(Pmat, vec)
!!    node_list%node(vg)%psi_eq(dof2) = vec(1)
!!    node_list%node(vg)%psi_eq(dof4) = vec(2)    
!!#endif
!
!    do in = i_tor_min, i_tor_max
!    do ivar = 1, n_var
!
!       vec(1) = node_list%node(vg)%values(in,dof2,ivar)
!       vec(2) = node_list%node(vg)%values(in,dof4,ivar)
!       vec    = matmul(Pmat, vec)
!    
!       node_list%node(vg)%values(in,dof2,ivar) = vec(1)
!       node_list%node(vg)%values(in,dof4,ivar) = vec(2)
!
!       vec(1) = node_list%node(vg)%deltas(in,dof2,ivar)
!       vec(2) = node_list%node(vg)%deltas(in,dof4,ivar)
!       vec    = matmul(Pmat, vec)
!
!       node_list%node(vg)%deltas(in,dof2,ivar) = vec(1)
!       node_list%node(vg)%deltas(in,dof4,ivar) = vec(2)
!
!    enddo
!    enddo
!
!  endif
!
!enddo
!end subroutine transform_nodelist
!
!subroutine transform_back_nodelist(node_list, i_tor_min, i_tor_max)
!use data_structure
!use phys_module
!implicit none
!integer :: vg, dof2, dof4, in, ivar
!type(type_node_list),    intent(inout) :: node_list
!integer ,              intent(in)    ::  i_tor_min, i_tor_max
!real*8  :: Pmat(2,2), vec(2), aa, bb, cc, dd
!
!dof2 = 2 ; dof4 = 4
!
!do vg = 1, node_list%n_nodes
!
!  if ( node_list%node(vg)%axis_node ) then
!
!    aa = node_list%node(vg)%x(1,dof2,1)  ;  bb = node_list%node(vg)%x(1,dof2,2)
!    cc = node_list%node(vg)%x(1,dof4,1)  ;  dd = node_list%node(vg)%x(1,dof4,2)
!          
!    Pmat(1,1) =  dd  ;  Pmat(1,2) = -bb
!    Pmat(2,1) = -cc  ;  Pmat(2,2) =  aa
!
!    Pmat = Pmat / (aa*dd - bb*cc)
!!#if fullmhd
!!    vec(1) = node_list%node(vg)%Fprof_eq(dof2)
!!    vec(2) = node_list%node(vg)%Fprof_eq(dof4)
!!    vec    = matmul(Pmat, vec)
!!    node_list%node(vg)%Fprof_eq(dof2) = vec(1)
!!    node_list%node(vg)%Fprof_eq(dof4) = vec(2) 
!!
!!    vec(1) = node_list%node(vg)%psi_eq(dof2)
!!    vec(2) = node_list%node(vg)%psi_eq(dof4)
!!    vec    = matmul(Pmat, vec)
!!    node_list%node(vg)%psi_eq(dof2) = vec(1)
!!    node_list%node(vg)%psi_eq(dof4) = vec(2)    
!!#elif altcs
!!    vec(1) = node_list%node(vg)%psi_eq(dof2)
!!    vec(2) = node_list%node(vg)%psi_eq(dof4)
!!    vec    = matmul(Pmat, vec)
!!    node_list%node(vg)%psi_eq(dof2) = vec(1)
!!    node_list%node(vg)%psi_eq(dof4) = vec(2)
!!#endif
!
!    do in = i_tor_min, i_tor_max
!    do ivar = 1, n_var
!
!       vec(1) = node_list%node(vg)%values(in,dof2,ivar)
!       vec(2) = node_list%node(vg)%values(in,dof4,ivar)
!       vec    = matmul(Pmat, vec)
!    
!       node_list%node(vg)%values(in,dof2,ivar) = vec(1)
!       node_list%node(vg)%values(in,dof4,ivar) = vec(2)
!
!       vec(1) = node_list%node(vg)%deltas(in,dof2,ivar)
!       vec(2) = node_list%node(vg)%deltas(in,dof4,ivar)
!       vec    = matmul(Pmat, vec)
!
!       node_list%node(vg)%deltas(in,dof2,ivar) = vec(1)
!       node_list%node(vg)%deltas(in,dof4,ivar) = vec(2)
!
!    enddo
!    enddo
!
!  endif
!
!enddo
!end subroutine transform_back_nodelist

end module mod_axis_treatment

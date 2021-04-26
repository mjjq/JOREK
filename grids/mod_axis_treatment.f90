module mod_axis_treatment

contains

subroutine new2old_dofs_on_the_axis(node_list, element, nodes)
use data_structure
implicit none
type(type_node_list),  intent(in)    :: node_list
type(type_element),    intent(inout) :: element
type(type_node),       intent(inout) :: nodes(n_vertex_max)
! --- routine parameters
integer :: iv, vg, ivar, in, dof2, dof4

dof2 = 2 ; dof4 = 4

do iv = 1, 4, 3  ! loop over axis nodes only, 1 and 4 are axis nodes
vg = element%vertex(iv)
  do in=1,n_tor
    do ivar = 1, n_var
      nodes(iv)%values(in,dof2,ivar) = nodes(iv)%x(1,dof2,1) * node_list%node(vg)%values(in,dof2,ivar) &
                                     + nodes(iv)%x(1,dof2,2) * node_list%node(vg)%values(in,dof4,ivar)
      nodes(iv)%values(in,dof4,ivar) = nodes(iv)%x(1,dof4,1) * node_list%node(vg)%values(in,dof2,ivar) &
                                     + nodes(iv)%x(1,dof4,2) * node_list%node(vg)%values(in,dof4,ivar)

      nodes(iv)%deltas(in,dof2,ivar) = nodes(iv)%x(1,dof2,1) * node_list%node(vg)%deltas(in,dof2,ivar) &
                                     + nodes(iv)%x(1,dof2,2) * node_list%node(vg)%deltas(in,dof4,ivar)
      nodes(iv)%deltas(in,dof4,ivar) = nodes(iv)%x(1,dof4,1) * node_list%node(vg)%deltas(in,dof2,ivar) &
                                     + nodes(iv)%x(1,dof4,2) * node_list%node(vg)%deltas(in,dof4,ivar)
    enddo
  enddo
enddo
end subroutine new2old_dofs_on_the_axis

subroutine old2new_basis_on_the_axis(nodes, element, ELM, RHS, i_tor_min, i_tor_max)
use data_structure
use phys_module, only: n_tor_fft_thresh
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
type(type_element),    intent(inout) :: element
#define DIM0 n_tor*n_vertex_max*(n_order+1)*n_var
real*8,          intent(inout) :: ELM(1:DIM0, 1:DIM0)
real*8,          intent(inout) :: RHS(1:DIM0)
integer ,        intent(in)    ::  i_tor_min, i_tor_max
! --- routine parameters
integer :: axis_vertex1, axis_vertex4, dof1, dof2, dof3, dof4
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo, ivar, jvar
integer :: im, in, n_tor_start, n_tor_end, n_tor_local
real*8  :: Ptrans(1:4, 1:4), Pmat(1:4, 1:4), rhs_i(1:4), elm_ij(1:4, 1:4)
logical :: use_fft

axis_vertex1 = 1 ; axis_vertex4 = 4
dof1 = 1 ; dof2 = 2 ; dof3 = 3 ; dof4 = 4

if ( (i_tor_min == 1) .and. (i_tor_max == n_tor) ) then
  ! In case of global matrix construction:
  use_fft = n_tor > n_tor_fft_thresh
else
  ! In case of "direct construction" of harmonic matrix never FFT:
  use_fft = .false.
end if
if ( use_fft ) then
  ! In case of FFT, don't loop over toroidal harmonics:
  n_tor_start = 1
  n_tor_end   = 1
else
  n_tor_start = i_tor_min
  n_tor_end   = i_tor_max
end if
n_tor_local = n_tor_end - n_tor_start +1

do iv = 1, n_vertex_max

  ! initialize traspose of the transformation matrix to identity
  Ptrans = 0.d0
  do io = 1, n_order+1
    Ptrans(io, io) = 1.d0
  enddo

  ! determine traspose of the transformation matrix
  if(iv==axis_vertex1 .or. iv==axis_vertex4)then
    Ptrans(dof2,dof2) = nodes(iv)%x(1,dof2,1) ; Ptrans(dof2,dof4) = nodes(iv)%x(1,dof4,1)
    PTrans(dof3,dof3) = 0.d0
    Ptrans(dof4,dof2) = nodes(iv)%x(1,dof2,2) ; Ptrans(dof4,dof4) = nodes(iv)%x(1,dof4,2)
  endif

  do im = n_tor_start, n_tor_end
  do ivar = 1, n_var  

    ! Extract RHS associated with a vertex iv: rhs_iv
    do io = 1, n_order+1
      if (use_fft) then
        index_iv_io =       n_var*(n_order+1)*(iv-1) +       n_var*(io-1) + 1 + (ivar-1)
      else
        index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + im - n_tor_start +1 + (ivar-1)*n_tor_local
      endif
      rhs_i(io)   = RHS(index_iv_io)
    enddo

    ! transform test functions for rhs_v
    rhs_i = matmul( Ptrans, rhs_i)

    ! fill the updated entries in RHS vector
    do io = 1, n_order+1
      if (use_fft) then
        index_iv_io =       n_var*(n_order+1)*(iv-1) +       n_var*(io-1) + 1 + (ivar-1)
      else
        index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + im - n_tor_start +1 + (ivar-1)*n_tor_local
      endif
      RHS(index_iv_io) = rhs_i(io)
    enddo

  enddo
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
    do im = n_tor_start, n_tor_end
    do ivar = 1, n_var
      do io = 1, n_order+1
        if (use_fft) then
          index_iv_io =       n_var*(n_order+1)*(iv-1) +       n_var*(io-1) + 1 + (ivar-1)
        else
          index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + im - n_tor_start +1 + (ivar-1)*n_tor_local
        endif
        do in =  n_tor_start, n_tor_end
        do jvar = 1, n_var
           do jo = 1, n_order+1
             if (use_fft) then
               index_jv_jo =       n_var*(n_order+1)*(jv-1) +       n_var*(jo-1) + 1 + (jvar-1)
             else
               index_jv_jo = n_tor_local*n_var*(n_order+1)*(jv-1) + n_tor_local*n_var*(jo-1) + in - n_tor_start +1 + (jvar-1) * n_tor_local
             endif
             elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
           enddo
        enddo
        enddo
      enddo
    enddo
    enddo

    ! transform test and basis functions for elm_iv_jv   
    elm_ij = matmul(matmul(PTrans, elm_ij), Pmat)

    ! fill the updated entries in ELM matrix
    do im = n_tor_start, n_tor_end
    do ivar = 1, n_var
      do io = 1, n_order+1
        if (use_fft) then
          index_iv_io =       n_var*(n_order+1)*(iv-1) +       n_var*(io-1) + 1 + (ivar-1)
        else
          index_iv_io = n_tor_local*n_var*(n_order+1)*(iv-1) + n_tor_local*n_var*(io-1) + im - n_tor_start +1 + (ivar-1)*n_tor_local
        endif
        do in =  n_tor_start, n_tor_end
        do jvar = 1, n_var
          do jo = 1, n_order+1
             if (use_fft) then
               index_jv_jo =       n_var*(n_order+1)*(jv-1) +       n_var*(jo-1) + 1 + (jvar-1)
             else
               index_jv_jo = n_tor_local*n_var*(n_order+1)*(jv-1) + n_tor_local*n_var*(jo-1) + in - n_tor_start +1 + (jvar-1) * n_tor_local
             endif
            ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
          enddo
        enddo
        enddo
      enddo
    enddo
    enddo

  enddo

enddo

end subroutine old2new_basis_on_the_axis

subroutine identify_axis_elements(node_list,element_list)
use data_structure
use phys_module
implicit none
integer :: ie, iv, j1, j2 ,j3 ,j4
type(type_node_list),    intent(inout) :: node_list
type(type_element_list), intent(inout) :: element_list

do ie = 1, element_list%n_elements
  j1 = element_list%element(ie)%vertex(1)
  j2 = element_list%element(ie)%vertex(2)
  j3 = element_list%element(ie)%vertex(3)
  j4 = element_list%element(ie)%vertex(4)

  element_list%element(ie)%axis_element = .false.

  ! The first and fourth vertex is on the grid-axis.
  if ( (treat_axis .or. treat_axis2) .and. ( node_list%node(j1)%axis_node .and. node_list%node(j4)%axis_node) ) then
     element_list%element(ie)%axis_element = .true.
  endif
enddo

end subroutine identify_axis_elements

end module mod_axis_treatment

module mod_assembly
implicit none
!> This module provides routines to factorize code in boundary conditions from each models

contains
  !>
  !! Subroutine: boundary_conditions_add_one_entry
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param solve_only  Do not add to harmonic matrix if .true.
  !! @param gmres       Do not add to product matrix if .false.
  !! @param cnt         Entry counter for precond murge problem.
  !! @param cnt_prod    Entry counter for product matrix.
  !! @param only_count  Indicate if we do a count or a real assembly.
  !! @param use_murge   Use murge interface.
  !! @param use_murge_element Murge interface with elementary matrices.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_one_entry( &
       &   index_node,  k,  in,                 &
       &   index_node2, k2, in2,                &
       &   zbig, solve_only, gmres,             &
       &   cnt, cnt_prod, only_count,           &
       &   index_min, index_max,                & 
       &   ijA_index, ijA_size, irn_jcn,        & 
       &   irn_glob, jcn_glob, A_glob,          & 
       &   i_tor_min, i_tor_max)
    use mod_parameters
!    use global_distributed_matrix
    use mod_locate_irn_jcn

    integer, intent(in)    :: index_node,  k,  in
    integer, intent(in)    :: index_node2, k2, in2
    real*8,  intent(in)    :: zbig
    logical, intent(in)    :: solve_only, gmres
    integer, intent(inout) :: cnt, cnt_prod
    logical, intent(in)    :: only_count
    integer, intent(in)    :: index_min, index_max
    integer, intent(in)    :: i_tor_min, i_tor_max 
    integer, intent(in), pointer :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer :: irn_glob(:), jcn_glob(:) 
    real*8  :: A_glob(:) 
    logical :: is_local
    integer :: ija_position, ilarge_vp

    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

       call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position,& 
                                    ijA_index, ijA_size, irn_jcn)
                             
       !-------- index dans A_glob
       ilarge_vp  = ijA_position  - 1 + ((k-1)*(i_tor_max - i_tor_min +1) + in-i_tor_min ) * n_var*(i_tor_max - i_tor_min +1) + (k2-1)*(i_tor_max - i_tor_min +1) + in2&
                    -i_tor_min + 1 
                               
                             
       irn_glob(ilarge_vp) =  (i_tor_max - i_tor_min +1) * n_var * (index_node -1) + (k -1)*(i_tor_max - i_tor_min +1) + in - i_tor_min + 1
       jcn_glob(ilarge_vp) =  (i_tor_max - i_tor_min +1) * n_var * (index_node2-1) + (k2-1)*(i_tor_max - i_tor_min +1) + in2 - i_tor_min + 1
       A_glob(ilarge_vp)   = ZBIG
    endif
  end subroutine boundary_conditions_add_one_entry

  !>
  !! Subroutine: boundary_conditions_add_RHS
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param solve_only  Do not add to harmonic matrix if .true.
  !! @param gmres       Do not add to product matrix if .false.
  !! @param cnt         Entry counter for precond murge problem.
  !! @param cnt_prod    Entry counter for product matrix.
  !! @param only_count  Indicate if we do a count or a real assembly.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_RHS( &
       &   index_node, k, in,             &
       &   index_min, index_max,          &
       &   rhs_loc,  val,                 &
       &   i_tor_min, i_tor_max)
    use mod_parameters
    integer, intent(in)    :: index_node,  k,  in
    integer, intent(in)    :: index_min, index_max
    integer, intent(in)    :: i_tor_min, i_tor_max 
    real*8,  intent(in)    :: val
    real*8,  intent(inOUT) :: rhs_loc(*)

    logical :: is_local

    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
       RHS_loc((i_tor_max - i_tor_min +1)*n_var * (index_node-1) + (k-1)*(i_tor_max - i_tor_min +1) + in - i_tor_min + 1) = val
    endif
  end subroutine boundary_conditions_add_RHS
end module mod_assembly

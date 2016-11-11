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
  SUBROUTINE boundary_conditions_add_one_entry( &
       &   index_node,  k,  in,                 &
       &   index_node2, k2, in2,                &
       &   zbig, solve_only, gmres,             &
       &   use_murge, use_murge_element,        &
       &   cnt, cnt_prod, only_count,           &
       &   index_min, index_max)
    use mod_parameters
    USE global_distributed_matrix
    USE murge_module, ONLY : murge_add_one_entry, vertex_is_local
    USE mod_locate_irn_jcn

    INTEGER, INTENT(IN)    :: index_node,  k,  in
    INTEGER, INTENT(IN)    :: index_node2, k2, in2
    REAL*8,  INTENT(IN)    :: zbig
    LOGICAL, INTENT(IN)    :: solve_only, gmres, use_murge, use_murge_element
    INTEGER, INTENT(INOUT) :: cnt, cnt_prod
    LOGICAL, INTENT(IN)    :: only_count
    INTEGER, INTENT(IN)    :: index_min, index_max

    LOGICAL :: is_local
    integer :: ijA_position, ilarge_vp

    if (use_murge .and. use_murge_element) then
       call vertex_is_local(index_node, is_local)
       if (is_local) then
          call murge_add_one_entry( index_node,  k,  in,  &
               &                    index_node2, k2, in2, &
               &                    zbig, solve_only,  &
               &                    gmres,             &
               &                    cnt, cnt_prod,     &
               &                    only_count)
       endif
    else
       if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

          call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position)
                                
          !-------- index dans A_glob
          ilarge_vp  = ijA_position  - 1 + ((k-1)*n_tor + in-1) * n_var*n_tor + (k2-1)*n_tor + in2
                                  
                                
          irn_glob(ilarge_vp) =  n_tor * n_var * (index_node -1) + (k -1)*n_tor + in
          jcn_glob(ilarge_vp) =  n_tor * n_var * (index_node2-1) + (k2-1)*n_tor + in2
          A_glob(ilarge_vp)   = ZBIG
       endif
    endif
  END SUBROUTINE boundary_conditions_add_one_entry

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
  !! @param use_murge   Use murge interface.
  !! @param use_murge_element Murge interface with elementary matrices.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  SUBROUTINE boundary_conditions_add_RHS( &
       &   index_node, k, in,             &
       &   use_murge, use_murge_element,  &
       &   index_min, index_max,          &
       &   RHS_loc,  val)
    use mod_parameters
    use murge_module, only : murge_harmonic, vertex_is_local
    INTEGER, INTENT(IN)    :: index_node,  k,  in
    LOGICAL, INTENT(IN)    :: use_murge, use_murge_element
    INTEGER, INTENT(IN)    :: index_min, index_max
    REAL*8,  INTENT(IN)    :: val
    real*8,  INTENT(INOUT) :: rhs_loc(*)

    LOGICAL :: is_local

    ! -- with murge, local_elms is replicated on each harmonic. we don't want to enter 
    ! -- boundary conditions twice
    if (use_murge .and. use_murge_element .and. murge_harmonic .eq. 1) then
       call vertex_is_local(index_node, is_local)
       if (is_local) then
          RHS_loc(n_tor*n_var * (index_node-1) + (k-1)*n_tor + in) = val
       endif
    else
       if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
          RHS_loc(n_tor*n_var * (index_node-1) + (k-1)*n_tor + in) = val
       endif
    endif
  END SUBROUTINE boundary_conditions_add_RHS
end module mod_assembly

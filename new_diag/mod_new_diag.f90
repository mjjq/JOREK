!> This module makes the whole new_diag package available and adds some easy to use functionality.
module mod_new_diag
  
  
  
  
  
  use parameters
  use mod_position
  use mod_straight_field_line
  use mod_expression
  use mod_four_filter
  use mod_diag_output
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  integer, parameter :: HIGHFIELD_SIDE = 0
  integer, parameter :: LOWFIELD_SIDE  = 1
  integer, parameter :: BOTH_SIDES     = 2
  
  
  
  
  
  
  contains
  
  
  
  
  
  !> Init???
  
  
  
  
  
  !> Toroidally averaged expressions on the midplane.
  subroutine midplane_profiles(node_list, element_list, eq, units, expr_list, res1d, side, n_pts,  &
    ierr)
    
    ! --- Routine parameters
    type(type_node_list),           intent(in)    :: node_list
    type(type_element_list),        intent(in)    :: element_list
    type(t_equil_state),            intent(in)    :: eq
    integer,                        intent(in)    :: units !< Output in which units?
    type(t_expr_list),              intent(in)    :: expr_list
    real*8, allocatable,            intent(inout) :: res1d(:,:)
    integer,                        intent(in)    :: side
    integer,                        intent(in)    :: n_pts
    integer,                        intent(out)   :: ierr
    
    ! --- Local variables
    real*8               :: Rstart, Rend
    real*8, allocatable  :: result(:,:,:,:)
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    
    ierr = 0
    
    if ( side == LOWFIELD_SIDE ) then
      Rstart = eq%R_axis     + 1.d-3
      Rend   = eq%R_midpl(2) - 1.d-3
    else if ( side == HIGHFIELD_SIDE ) then
      Rstart = eq%R_midpl(1) + 1.d-3
      Rend   = eq%R_axis     - 1.d-3
    else if ( side == BOTH_SIDES ) then
      Rstart = eq%R_midpl(1) + 1.d-3
      Rend   = eq%R_midpl(2) - 1.d-3
    else
      !### error
    end if
    pol_pos_list = pol_pos(node_list, element_list, eq, Rstart=Rstart, Rend=Rend, Z=eq%Z_axis,     &
      n=n_pts)
    tor_pos_list = tor_pos(nphi=4*n_plane) !###
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call transform_and_filter(result, simple_filter(n=0))
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    deallocate(result)
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)
    
  end subroutine midplane_profiles
  
  
  
  
  
  !> Construct poloidally and toroidally averaged profiles.
  !subroutine average_profiles(eq, res1d, ierr)
  !end subroutine average_profiles
  
  
  
  
  
end module mod_new_diag

!> Module for general purpose diagnostic routines.
module diagnostics
  
  use mod_parameters
  use nodes_elements
  use phys_module
  use mod_interp
  
  
  
  implicit none
  
  
  
  type :: t_position_buffer
    real*8  :: R, Z, s, t
    integer :: ielm
  end type t_position_buffer
  
  
  
  !> This interface allows to determine the value of a variable at a given position
  !! easily. If many evaluations are necessary, it is not as performant as a direct
  !! reconstruction, but in many cases it should be sufficient.
  interface variable_value
    module procedure variable_value_R_Z_phi
    module procedure variable_value_ielm_s_t
    module procedure variable_value_inode
  end interface variable_value
  
  
  
  contains
  
  
  
  !> Determines the value of variable ivar at position (R,Z,phi)
  real*8 function variable_value_R_Z_phi(ivar, R, Z, phi, no_zero, error)
    
    ! --- Routine parameters
    integer, intent(in) :: ivar      !< Number of the variable
    real*8,  intent(in) :: R         !< R-position to evaluate the variable at
    real*8,  intent(in) :: Z         !< Z-position to evaluate the variable at
    real*8,  intent(in) :: phi       !< phi-position to evaluate the variable at   
    logical, intent(in) :: no_zero   !< Skip n=0 part?
    integer, intent(out):: error     !< Error flag

    ! --- Local variables
    logical,                 save :: initialized = .false.
    type(t_position_buffer), save :: position_buffer !< Buffer to avoid unnecessary find_RZ calls
    integer :: ielm, i_harm
    real*8  :: R_out, Z_out, s, t, P, P_s, P_t, P_st, P_ss, P_tt
    
    error = 0
    
    if ( .not. initialized ) then
      initialized       = .true.
      position_buffer%R = -99.d0
      position_buffer%Z = -99.d0
    end if
    
    if ( (R == position_buffer%R) .and. (Z == position_buffer%Z) ) then
      s    = position_buffer%s
      t    = position_buffer%t
      ielm = position_buffer%ielm
    else
      call find_RZ(node_list, element_list, R, Z, R_out, Z_out, ielm, s, t, error)
      if ( error /= 0 ) return
      position_buffer%R    = R
      position_buffer%Z    = Z
      position_buffer%s    = s
      position_buffer%t    = t
      position_buffer%ielm = ielm
    end if
    
    variable_value_R_Z_phi = variable_value_ielm_s_t(ivar, ielm, s, t, phi, no_zero)
    
  end function variable_value_R_Z_phi
  
  
  
  !> Determines the value of variable ivar at position (s,t,phi) in element i_elm
  !> This function can be replaced with interp_0 if no_zero == .false.
  real*8 function variable_value_ielm_s_t(ivar, i_elm, s, t, phi, no_zero)
    
    ! --- Routine parameters
    integer, intent(in) :: ivar      !< Number of the variable
    integer, intent(in) :: i_elm     !< Element index
    real*8,  intent(in) :: s         !< s-position to evaluate the variable at
    real*8,  intent(in) :: t         !< t-position to evaluate the variable at
    real*8,  intent(in) :: phi       !< phi-position to evaluate the variable at   
    logical, intent(in) :: no_zero   !< Skip n=0 part?
    
    ! --- Local variables
    integer :: i_harm
    real*8  :: P, P_s, P_t, P_st, P_ss, P_tt
    
    variable_value_ielm_s_t = 0.d0
    do i_harm = 1, n_tor
      
      if ( no_zero .and. (i_harm==1) ) cycle
      
      call interp(node_list, element_list, i_elm, ivar, i_harm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
      
      if ( i_harm == 1 ) then
         variable_value_ielm_s_t = variable_value_ielm_s_t + P
      else if ( mod(i_harm,2) == 0 ) then
         variable_value_ielm_s_t = variable_value_ielm_s_t + P * cos( mode(i_harm) * phi )
      else
         variable_value_ielm_s_t = variable_value_ielm_s_t + P * sin( mode(i_harm) * phi )
      end if
      
    end do
    
  end function variable_value_ielm_s_t
  
  
  
  !> Determines the value of variable ivar at toroidal position phi and node i_node
  real*8 function variable_value_inode(ivar, i_node, phi, no_zero)
    
    ! --- Routine parameters
    integer, intent(in) :: ivar      !< Number of the variable
    integer, intent(in) :: i_node    !< Node index
    real*8,  intent(in) :: phi       !< phi-position to evaluate the variable at   
    logical, intent(in) :: no_zero   !< Skip n=0 part?
    
    ! --- Local variables
    integer :: i_harm
    real*8  :: P, P_s, P_t, P_st, P_ss, P_tt
    
    variable_value_inode = 0.d0
    do i_harm = 1, n_tor
      
      if ( no_zero .and. (i_harm==1) ) cycle
      
      P = node_list%node(i_node)%values(i_harm, 1, ivar)
      
      if ( i_harm == 1 ) then
         variable_value_inode = variable_value_inode + P
      else if ( mod(i_harm,2) == 0 ) then
         variable_value_inode = variable_value_inode + P * cos( mode(i_harm) * phi )
      else
         variable_value_inode = variable_value_inode + P * sin( mode(i_harm) * phi )
      end if
      
    end do
    
  end function variable_value_inode
  
  
  
  !> Detects if the poloidal flux has a minimum or a maximum at the magnetic axis
  !! (i.e., the sign of the plasma current).
  logical function axis_is_psi_minimum(node_list, element_list, R_axis, Z_axis, psi_axis)
    
    use data_structure
    
    implicit none
    
    ! --- Routine parameters
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    real*8,                   intent(in) :: R_axis
    real*8,                   intent(in) :: Z_axis
    real*8,                   intent(in) :: Psi_axis
    
    ! --- Local variables
    integer :: i_elm, ifail, i_var, i_harm
    real*8  :: R_out, Z_out, s, t, P, P_s, P_t, P_st, P_ss, P_tt
    
    call find_RZ(node_list, element_list, R_axis+0.01, Z_axis, R_out, Z_out, i_elm, s, t, ifail)
    
    call interp(node_list, element_list, i_elm, 1, 1, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
    
    axis_is_psi_minimum = ( psi_axis < P )
    
  end function axis_is_psi_minimum
  
  
  
end module diagnostics

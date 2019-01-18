!> Module which allows to find out if a certain position is located in the plasma, the scrape-off
!! layer or a private flux region.
module domains
  
  
  
  use data_structure
  
  
  
  implicit none
  
  
  
  private
  public which_domain, in_plasma, in_sol, in_private
  
  
  
  
  logical, save :: axis_is_min
  logical, save :: axis_is_min_initialized = .false.
  
  
  
  contains
    
    
    
  !> Determine in which domain a certain position is located (plasma, vacuum, private fluxregion).
  integer function which_domain(node_list, element_list, R, Z, psi, xpoint, xcase, R_xpoint,       &
    Z_xpoint, psi_xpoint, psi_limit, R_axis, Z_axis, psi_axis)
    
    use phys_module, only: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL, DOMAIN_UPPER_PRIVATE,        &
      DOMAIN_LOWER_PRIVATE, LOWER_XPOINT, UPPER_XPOINT, DOUBLE_NULL
    use diagnostics, only: axis_is_psi_minimum
    
    type(type_node_list),    intent(in)  :: node_list    !< List of grid nodes
    type(type_element_list), intent(in)  :: element_list !< List of grid elements
    real*8,  intent(in) :: R             !< R position
    real*8,  intent(in) :: Z             !< Z position
    real*8,  intent(in) :: psi           !< Poloidal flux at (R,Z)
    logical, intent(in) :: xpoint        !< X-point case?
    integer, intent(in) :: xcase         !< Lower/upper X-point or double-null?
    real*8,  intent(in) :: R_xpoint(2)   !< R-position of X-point(s)
    real*8,  intent(in) :: Z_xpoint(2)   !< Z-position of X-point(s)
    real*8,  intent(in) :: psi_xpoint(2) !< Psi-value(s) at X-point(s)
    real*8,  intent(in) :: psi_limit     !< psi-value at limiter
    real*8,  intent(in) :: R_axis        !< R-position of mag. axis
    real*8,  intent(in) :: Z_axis        !< Z-position of mag. axis
    real*8,  intent(in) :: psi_axis      !< psi-value at axis
    
    real*8 :: sign_corr

    if ( .not. axis_is_min_initialized ) then
      axis_is_min = axis_is_psi_minimum(node_list, element_list, R_axis, Z_axis, psi_axis)
      axis_is_min_initialized = .true.
    end if
    
    if ( axis_is_min ) then
      sign_corr = 1.
    else
      sign_corr = -1.
    end if
 
    if ( xpoint ) then
      
      if ( xcase == LOWER_XPOINT ) then
        
        if ( ( sign_corr*psi < sign_corr*psi_xpoint(1) ) .and. ( Z < Z_xpoint(1) ) ) then
          which_domain = DOMAIN_LOWER_PRIVATE
        else if ( sign_corr*psi < sign_corr*psi_xpoint(1) ) then
          which_domain = DOMAIN_PLASMA
        else
          which_domain = DOMAIN_SOL
        end if
        
      else if ( xcase == UPPER_XPOINT ) then
        
        if ( ( sign_corr*psi < sign_corr*psi_xpoint(1) ) .and. ( Z > Z_xpoint(1) ) ) then
          which_domain = DOMAIN_UPPER_PRIVATE
        else if ( sign_corr*psi < sign_corr*psi_xpoint(1) ) then
          which_domain = DOMAIN_PLASMA
        else
          which_domain = DOMAIN_SOL
        end if
        
      else if ( xcase == DOUBLE_NULL ) then
        
        if ( ( sign_corr*psi < sign_corr*psi_xpoint(1) ) .and. ( Z < Z_xpoint(1) ) ) then
          which_domain = DOMAIN_LOWER_PRIVATE
        else if ( ( sign_corr*psi < sign_corr*psi_xpoint(2) ) .and. ( Z > Z_xpoint(2) ) ) then
          which_domain = DOMAIN_UPPER_PRIVATE
        else if ( sign_corr*psi < minval(sign_corr*psi_xpoint) ) then
          which_domain = DOMAIN_PLASMA
        else if ( sign_corr*psi < maxval(sign_corr*psi_xpoint) ) then
          which_domain = DOMAIN_SOL
        else
          which_domain = DOMAIN_OUTER_SOL
        end if
        
      end if
      
    else
      
      if ( sign_corr*psi < sign_corr*psi_limit ) then
        which_domain = DOMAIN_PLASMA
      else
        which_domain = DOMAIN_SOL
      end if
      
    end if
    
  end function which_domain
  

  !> Determine if a position is located inside the plasma
  logical function in_plasma(node_list, element_list, R, Z, psi, xpoint, xcase, R_xpoint, Z_xpoint,&
    psi_xpoint, psi_limit, R_axis, Z_axis, psi_axis)
    
    use phys_module, only: DOMAIN_PLASMA
    
    type(type_node_list),    intent(in)  :: node_list    !< List of grid nodes
    type(type_element_list), intent(in)  :: element_list !< List of grid elements
    real*8,  intent(in) :: R             !< R position
    real*8,  intent(in) :: Z             !< Z position
    real*8,  intent(in) :: psi           !< Poloidal flux at (R,Z)
    logical, intent(in) :: xpoint        !< X-point case?
    integer, intent(in) :: xcase         !< Lower/upper X-point or double-null?
    real*8,  intent(in) :: R_xpoint(2)   !< R-position of X-point(s)
    real*8,  intent(in) :: Z_xpoint(2)   !< Z-position of X-point(s)
    real*8,  intent(in) :: psi_xpoint(2) !< Psi-value(s) at X-point(s)
    real*8,  intent(in) :: psi_limit     !< psi-value at limiter
    real*8,  intent(in) :: R_axis        !< R-position of mag. axis
    real*8,  intent(in) :: Z_axis        !< Z-position of mag. axis
    real*8,  intent(in) :: psi_axis      !< psi-value at axis
    
    integer :: domain

    domain = which_domain(node_list,element_list,R,Z,psi,xpoint,xcase,R_xpoint,Z_xpoint,psi_xpoint,&
      psi_limit,R_axis,Z_axis,psi_axis)

    in_plasma = ( domain==DOMAIN_PLASMA )
 
  end function in_plasma
  
  
  
  !> Determine if a position is located in the scrape-off layer
  logical function in_sol(node_list, element_list, R, Z, psi, xpoint, xcase, R_xpoint, Z_xpoint,   &
    psi_xpoint, psi_limit, R_axis, Z_axis, psi_axis)
    
    use phys_module, only: DOMAIN_SOL, DOMAIN_OUTER_SOL
    
    type(type_node_list),    intent(in)  :: node_list    !< List of grid nodes
    type(type_element_list), intent(in)  :: element_list !< List of grid elements
    real*8,  intent(in) :: R             !< R position
    real*8,  intent(in) :: Z             !< Z position
    real*8,  intent(in) :: psi           !< Poloidal flux at (R,Z)
    logical, intent(in) :: xpoint        !< X-point case?
    integer, intent(in) :: xcase         !< Lower/upper X-point or double-null?
    real*8,  intent(in) :: R_xpoint(2)   !< R-position of X-point(s)
    real*8,  intent(in) :: Z_xpoint(2)   !< Z-position of X-point(s)
    real*8,  intent(in) :: psi_xpoint(2) !< Psi-value(s) at X-point(s)
    real*8,  intent(in) :: psi_limit     !< psi-value at limiter
    real*8,  intent(in) :: R_axis        !< R-position of mag. axis
    real*8,  intent(in) :: Z_axis        !< Z-position of mag. axis
    real*8,  intent(in) :: psi_axis      !< psi-value at axis
    
    integer :: domain
    
    domain = which_domain(node_list,element_list,R,Z,psi,xpoint,xcase,R_xpoint,Z_xpoint,psi_xpoint,&
      psi_limit,R_axis,Z_axis,psi_axis)
    
    in_sol = ( (domain==DOMAIN_SOL) .or. (domain==DOMAIN_OUTER_SOL) )
  
  end function in_sol
  
  
  
  !> Determine if a position is located in a private flux region
  logical function in_private(node_list, element_list, R, Z, psi, xpoint, xcase, R_xpoint,         &
    Z_xpoint, psi_xpoint, psi_limit, R_axis, Z_axis, psi_axis)
    
    use phys_module, only: DOMAIN_UPPER_PRIVATE, DOMAIN_LOWER_PRIVATE
    
    type(type_node_list),    intent(in)  :: node_list    !< List of grid nodes
    type(type_element_list), intent(in)  :: element_list !< List of grid elements
    real*8,  intent(in) :: R             !< R position
    real*8,  intent(in) :: Z             !< Z position
    real*8,  intent(in) :: psi           !< Poloidal flux at (R,Z)
    logical, intent(in) :: xpoint        !< X-point case?
    integer, intent(in) :: xcase         !< Lower/upper X-point or double-null?
    real*8,  intent(in) :: R_xpoint(2)   !< R-position of X-point(s)
    real*8,  intent(in) :: Z_xpoint(2)   !< Z-position of X-point(s)
    real*8,  intent(in) :: psi_xpoint(2) !< Psi-value(s) at X-point(s)
    real*8,  intent(in) :: psi_limit     !< psi-value at limiter
    real*8,  intent(in) :: R_axis        !< R-position of mag. axis
    real*8,  intent(in) :: Z_axis        !< Z-position of mag. axis
    real*8,  intent(in) :: psi_axis      !< psi-value at axis
    
    integer :: domain
    
    domain = which_domain(node_list,element_list,R,Z,psi,xpoint,xcase,R_xpoint,Z_xpoint,psi_xpoint,&
      psi_limit,R_axis,Z_axis,psi_axis)
    
    in_private = ( (domain==DOMAIN_UPPER_PRIVATE) .or. (domain==DOMAIN_LOWER_PRIVATE) )
  
  end function in_private
  
  
  
end module domains

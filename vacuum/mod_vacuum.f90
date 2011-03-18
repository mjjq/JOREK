!> Common parameters and variables for the JOREK free boundary extension.
!!
!! @see vacuum_response, vacuum_equilibrium
module vacuum
  
  
  implicit none
  
  
  !> @name General parameters
  logical, parameter  :: vacuum_debug          = .true.  !< Enable additional output and tests
  logical, parameter  :: vacuum_implicit       = .true.  !< Enable implicit time-evolution of wall currents @todo REMOVE
  logical, parameter  :: vacuum_decouple_modes = .false. !< Option to switch off 3D wall mode coupling
  logical, parameter  :: NEW_VACUUM            = .true.  !< @todo REMOVE
  integer             :: n_dof_bnd                       !< Total number of boundary dofs per harmonic
  integer             :: n_dof_starwall                  !< Total number of boundary dofs in STARWALL response
  
  !> @name Resistive wall only
  real*8              :: wall_resistivity    = 0.d0      !< Resistivity of the external wall @todo Units etc.???
  real*8              :: wall_thickness      = 1.d0      !< Thickness of the external wall @todo Units etc.???
  integer             :: n_wall_curr                     !< Number of wall current potentials.
  real*8, allocatable :: wall_curr(:)                    !< Wall current potentials (\f$Y_k\f$ in the documentation).

  !> @name STARWALL vacuum response
  integer              :: n_starwall_harmonics           !< Number of toroidal harmonics in response
  integer, allocatable :: starwall_harmonics(:)          !< Harmonics in the response
  real*8,  allocatable :: starwall_d_yy(:)               !< YY-matrix, see documentation
  real*8,  allocatable :: starwall_m_ye(:,:)             !< YE-matrix, see documentation
  real*8,  allocatable :: starwall_m_ey(:,:)             !< EY-matrix, see documentation
  real*8,  allocatable :: starwall_m_ee(:,:)             !< EE-matrix, see documentation
  real*8,  allocatable :: starwall_m_id(:,:)             !< Ideal wall matrix, see documentation
  real*8,  allocatable :: starwall_m_nw(:,:)             !< No wall matrix, see documentation

  !> @name JOREK vacuum response matrices
  !! Response matrices derived from STARWALL response (w=wall, p=plasma)
  real*8, allocatable :: response_d_ww(:)                !< \f$R_{k k}\f$ in the documentation
  real*8, allocatable :: response_m_wp(:,:)              !< \f$S_{k j}\f$ in the documentation
  real*8, allocatable :: response_m_pp_lhs(:,:)          !< \f$T_{i j}\f$ in the documentation
  real*8, allocatable :: response_m_pp_rhs(:,:)          !< \f$MMM_{i j}\f$ in the documentation
  real*8, allocatable :: response_m_pw_rhs(:,:)          !< \f$U_{i k}\f$ in the documentation
  real*8, allocatable :: response_m_eq(:,:)              !< Response matrix for vacuum_equil

 !> @name Coil contributions
  integer             :: n_coils                         !< number of poloidal field coils
  real*8, allocatable :: R_coils(:), Z_coils(:)          !< positions of poloidal field coils
  real*8, allocatable :: dR_coils(:), dZ_coils(:)        !< width/height of poloidal field coils
  real*8, allocatable :: I_coils(:)                      !< coil currents 
  real*8, allocatable :: external_field(:,:)             !< external poloidal field (n_dof_bnd,n_coils)
  
  !> @name OLD PARAMETERS -- WILL BE REMOVED SOON
  character(len=8)    :: wall_curr_treatment = 'implicit'!< @todo REMOVE (OLD)
  real*8, allocatable :: vac_response(:,:)               !< @todo REMOVE (OLD)
  real*8, allocatable :: vac_response2(:,:,:)            !< @todo REMOVE (OLD)
  real*8, allocatable :: vac_response3(:,:,:)            !< @todo REMOVE (OLD)
  real*8, allocatable :: diagonal_yy(:)                  !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_ye(:,:)                  !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_ey(:,:)                  !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_ee(:,:)                  !< @todo REMOVE (OLD)
  real*8, allocatable :: diagonal_r(:)                   !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_s(:,:)                   !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_t(:,:)                   !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_u(:,:)                   !< @todo REMOVE (OLD)
  real*8, allocatable :: matrix_v(:,:)                   !< @todo REMOVE (OLD)
  
  
  
  contains
  
  
  
  !> Preset freeboundary related input parameters to reasonable default values.
  subroutine vacuum_preset(my_id, freeboundary_equil, freeboundary, use_starwall, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(out)   :: freeboundary_equil, freeboundary, use_starwall, resistive_wall
    
    if ( my_id /= 0 ) return
    
    freeboundary_equil = .false. ! Free or fixed boundary equilibrium?
    freeboundary       = .false. ! Free or fixed boundary conditions in time-evolution?
    use_starwall       = .false. ! Use STARWALL vacuum solution?
    resistive_wall     = .false. ! Resistive or ideal wall?
    
  end subroutine vacuum_preset
  
  
  
  !> Initialize vacuum, ensure consistency of input parameters.
  subroutine vacuum_init(my_id, freeboundary_equil, freeboundary, use_starwall, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(inout) :: freeboundary_equil, freeboundary, use_starwall, resistive_wall
    
    if ( my_id /= 0 ) return
    
    ! --- Make input parameters consistent.
    if ( freeboundary_equil ) freeboundary = .true.
    if ( freeboundary) then
      if ( resistive_wall ) use_starwall = .true.
    else
      use_starwall   = .false.
      resistive_wall = .false.
    end if
        
  end subroutine vacuum_init

end module vacuum

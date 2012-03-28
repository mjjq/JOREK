!> Common parameters and variables for the JOREK free boundary extension.
!!
!! @see vacuum_response, vacuum_equilibrium
module vacuum
  
  implicit none
  
  !> @name General parameters
  logical, parameter  :: vacuum_debug          = .false. !< Enable additional output and tests
  logical, parameter  :: vacuum_decouple_modes = .false. !< Option to switch off 3D wall mode coupling
  integer             :: n_dof_bnd                       !< Total number of boundary dofs per harmonic
  integer             :: n_dof_starwall                  !< Total number of boundary dofs in STARWALL response
  integer, parameter  :: ivar_psi = 1                    !< Index of Psi variable
  integer, parameter  :: ivar_j   = 3                    !< Index of j variable
  
  !> @name Resistive wall only
  real*8              :: wall_resistivity    = 1.d-2     !< Resistivity of the external wall @todo Units etc.??? @todo Make input parameter
  real*8              :: wall_thickness      = 1.d0      !< Thickness of the external wall @todo Units etc.??? @todo Make input parameter
  integer             :: n_wall_curr                     !< Number of wall current potentials.
  real*8, allocatable :: wall_curr(:)                    !< Wall current potentials (\f$Y_k\f$).
  real*8, allocatable :: dwall_curr(:)                   !< Change of wall current potentials (\f$\delta Y_k\f$).
  real*8, allocatable :: old_dpsibnd_vec(:)              !< Previous delta Psi values required for time-stepping with zeta/=0

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
  real*8, allocatable :: response_m_a(:,:)               !< \f$\hat{A}\f$ in the documentation
  real*8, allocatable :: response_d_b(:)                 !< \f$\hat{B}\f$ in the documentation
  real*8, allocatable :: response_d_c(:)                 !< \f$\hat{C}\f$ in the documentation
  real*8, allocatable :: response_m_d(:,:)               !< \f$\hat{D}\f$ in the documentation
  real*8, allocatable :: response_m_e(:,:)               !< \f$\hat{E}\f$ in the documentation
  real*8, allocatable :: response_m_f(:,:)               !< \f$\hat{F}\f$ in the documentation
  real*8, allocatable :: response_m_g(:,:)               !< \f$\hat{G}\f$ in the documentation
  real*8, allocatable :: response_m_h(:,:)               !< \f$\hat{H}\f$ in the documentation
  real*8, allocatable :: response_m_j(:,:)               !< \f$\hat{J}\f$ in the documentation
  real*8, allocatable :: response_m_eq(:,:)              !< Response matrix for vacuum_equil

 !> @name Coil contributions
  integer             :: n_coils                         !< number of poloidal field coils
  real*8, allocatable :: R_coils(:), Z_coils(:)          !< positions of poloidal field coils
  real*8, allocatable :: dR_coils(:), dZ_coils(:)        !< width/height of poloidal field coils
  real*8, allocatable :: I_coils(:)                      !< coil currents 
  real*8              :: vertical_FB                     !< a variable for the feedback control of the plasma's vertical position
  real*8, allocatable :: bext_par(:,:)                   !< external poloidal field tangential to
                                                         !! the interface (n_dof_bnd,n_coils)
  
  
  
  contains
  
  
  
  !> Preset freeboundary related input parameters to reasonable default values.
  subroutine vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(out)   :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Preset namelist input parameters.
    freeboundary_equil   = .false.
    freeboundary         = .false.
    resistive_wall       = .false.
    
  end subroutine vacuum_preset
  
  
  
  !> Initialize vacuum, ensure consistency of input parameters.
  subroutine vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(inout) :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Make input parameters consistent.
    freeboundary   = freeboundary .or. freeboundary_equil
    resistive_wall = freeboundary .and. resistive_wall
        
    ! --- Initialize some variables.
    n_starwall_harmonics = 0
    n_wall_curr          = 0
    n_dof_bnd            = 0
    n_dof_starwall       = 0
    n_coils              = 0
    
  end subroutine vacuum_init
  
  
  
  !> Allows to decide if free- or fixed-boundary conditions apply for a certain combination of
  !! toroidal harmonic i_tor and variable i_var. This is important, as the STARWALL response
  !! must not be provided for all toroidal modes and fixed boundary conditions remain implemented
  !! for some of the JOREK variables.
  logical function is_freebound(i_tor, i_var)
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor
    integer, intent(in) :: i_var
    
    ! --- Local variables
    integer :: i
    
    ! --- Free boundary conditions only if STARWALL response provided for this toroidal harmonic
    is_freebound = .false.
    do i = 1, n_starwall_harmonics
      is_freebound = is_freebound .or. ( i_tor == starwall_harmonics(i) )
    end do
    
    ! --- Free boundary conditions only for certain variables
    is_freebound = is_freebound .and. ( (i_var == ivar_j) .or. (i_var == ivar_psi) )
    
  end function is_freebound
  
  
  
  !> Import some vacuum-related data from the restart file  
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_restart_vacuum(file_handle, freeboundary, resistive_wall, index_start)
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    integer, intent(in) :: index_start
    
    ! --- Local variables
    logical :: resistive_wall_rst
    
    if ( freeboundary .and. (index_start > 0) ) then
      
      read(21) resistive_wall_rst
      if ( resistive_wall /= resistive_wall_rst ) then
        write(*,*) 'ERROR: It is currently not possible to restart a JOREK simulation with a'
        write(*,*) '  modified setting for resistive_wall.'
        stop
      end if
      
      if ( resistive_wall ) then
        
        if ( allocated(wall_curr) ) deallocate(wall_curr)
        allocate( wall_curr(n_wall_curr) )
        read(21) wall_curr(:)
        
        if ( allocated(dwall_curr) ) deallocate(dwall_curr)
        allocate( dwall_curr(n_wall_curr) )
        read(21) dwall_curr(:)
        
      end if
      
    end if
    
    if ( vacuum_debug ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine import_restart_vacuum
  
  
  
  !> Export some vacuum-related data to the restart file  
  subroutine export_restart_vacuum(file_handle, freeboundary, resistive_wall, index_now)
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    integer, intent(in) :: index_now
    
    if ( freeboundary .and. (index_now > 0) ) then
      
      write(21) resistive_wall
      if ( resistive_wall ) then
        
        if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or.                    &
          (.not. allocated(old_dpsibnd_vec)) ) then
          write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
          stop
        end if
        
        write(file_handle) wall_curr(:)
        
        write(file_handle) dwall_curr(:)
        
      end if
      
    end if
    
    if ( vacuum_debug ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine export_restart_vacuum
  
  
  
end module vacuum

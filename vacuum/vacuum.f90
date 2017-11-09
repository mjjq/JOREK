!> Common parameters and variables for the JOREK free boundary extension.
!!
!! @see vacuum_response, vacuum_equilibrium
module vacuum
  
  implicit none
  
  !> @name General parameters
  logical, parameter  :: vacuum_debug          = .true. !< Enable additional output and tests
  logical, parameter  :: vacuum_decouple_modes = .false. !< Option to switch off 3D wall mode coupling
  integer             :: n_dof_bnd                       !< Total number of boundary dofs per harmonic
  integer             :: n_dof_starwall                  !< Total number of boundary dofs in STARWALL response
  integer, parameter  :: ivar_psi = 1                    !< Index of Psi variable
  integer, parameter  :: ivar_j   = 3                    !< Index of j variable
  real*8              :: freeb_fact                      !< Switches on free-boundary terms in elt_matrix when =1.
  
  !> @name Resistive wall only
  real*8              :: wall_resistivity_fact           !< Scaling factor for the wall resistivity specified in STARWALL
  real*8              :: wall_resistivity                !< Resistivity of the external wall
  logical             :: wall_curr_initialized = .false. !< Have the wall currents been initialized?
  integer             :: n_wall_curr                     !< Number of wall current potentials.
  real*8, allocatable :: wall_curr(:)                    !< Wall current potentials (\f$Y_k\f$).
  real*8, allocatable :: dwall_curr(:)                   !< Change of wall current potentials (\f$\delta Y_k\f$).
  real*8, allocatable :: old_dpsibnd_vec(:)              !< Previous delta Psi values required for time-stepping with zeta/=0

  !> Data type for response matrices distributed over the MPI tasks
  type :: t_distrib_mat
    real*8, allocatable :: loc_mat(:,:)                  !< Local chunk of the matrix
    logical             :: distrib                       !< Is the matrix distributed?
    logical             :: row_wise                      !< Is the matrix distributed rowwise (otherwise columnwise)?
    integer             :: ind_start                     !< Minimum row/column index of local chunk.
    integer             :: ind_end                       !< Maximum row/column index of local chunk.
    integer             :: step                          !< 
  end type t_distrib_mat

  !> @name JOREK vacuum response matrices
  !! Response matrices derived from STARWALL response (w=wall, p=plasma)
  type(t_distrib_mat)  :: response_m_a                   !< \f$\hat{A}\f$ in the documentation
  real*8, allocatable  :: response_d_b(:)                !< \f$\hat{B}\f$ in the documentation
  real*8, allocatable  :: response_d_c(:)                !< \f$\hat{C}\f$ in the documentation
  type(t_distrib_mat)  :: response_m_d                   !< \f$\hat{D}\f$ in the documentation
  real*8, allocatable  :: response_m_e(:,:)              !< \f$\hat{E}\f$ in the documentation
  type(t_distrib_mat)  :: response_m_f                   !< \f$\hat{F}\f$ in the documentation
  type(t_distrib_mat)  :: response_m_g                   !< \f$\hat{G}\f$ in the documentation
  real*8, allocatable  :: response_m_h(:,:)              !< \f$\hat{H}\f$ in the documentation
  real*8, allocatable  :: response_m_j(:,:)              !< \f$\hat{J}\f$ in the documentation
  real*8, allocatable  :: response_m_k(:,:)              !< \f$\hat{K}\f$ in the documentation
  real*8, allocatable  :: response_m_l(:,:)              !< \f$\hat{L}\f$ in the documentation
  type(t_distrib_mat)  :: response_m_v                   !< \f$\hat{V}\f$ in the documentation
  real*8, allocatable  :: response_m_eq(:,:)             !< Response matrix for vacuum_equil

  !> @name Equilibrium coil contributions
  integer             :: n_coils                         !< number of poloidal field coils in coil_field.dat
  integer             :: n_pf_coils                      !< number of poloidal field coils
  logical             :: starwall_equil_coils            !< specify wheter the equilibrium PF coils will be given by STARWALL or not
  real*8, allocatable :: I_coils(:)                      !< coil currents
  real*8, allocatable :: Y_coils0(:)                     !< imposed STARWALL coil currents source
  real*8              :: vertical_FB                     !< a variable for the feedback control of the plasma's vertical position
  real*8, allocatable :: bext_tan(:,:)                   !< external tangential field
  real*8, allocatable :: bext_nor(:,:)                   !< external normal field
  real*8, allocatable :: bext_psi(:,:)                   !< external poloidal flux      
  real*8              :: psi_offset_freeb                !< Allows to shift the value of psi by a global constant for freeb_equil (improves convergence)
  
  !> @name Equilibrium parameters for feedback
  real*8              :: current_ref                     !< Target total plasma current Ip for the feedback (FB)
  real*8              :: FB_Ip_position                  !< Amplification factor for Ip feedback
  real*8              :: FB_Ip_integral                  !< Amplification factor for Ip feedback
  real*8              :: Z_axis_ref                      !< Target magnetic axis vertical position
  real*8              :: FB_Zaxis_position               !< Amplification factor for Zaxis feedback
  real*8              :: FB_Zaxis_derivative             !< Amplification factor for Zaxis feedback
  real*8              :: FB_Zaxis_integral               !< Amplification factor for Zaxis feedback
  integer             :: start_VFB                       !< Iteration for starting vertical feedback
  integer             :: n_feedback_current              !< Feedback will be performed each n_... iterations
  integer             :: n_feedback_vertical             !< Feedback will be performed each n_... iterations
  integer             :: n_iter_freeb                    !< Number of iterations for freeboundary equilibirum
  
  !> @name Time-evolution PF coils parameters
  real*8              :: PF_pert_start_time              !< Time to start a perturbation to speed-up VDEs
  
  
  ! ### various variables, some need to be removed
  real*8, allocatable :: R_coils(:), Z_coils(:)          ! ### old
  real*8, allocatable :: dR_coils(:), dZ_coils(:)        ! ### old
  real*8, allocatable :: coil_voltages(:)                !< Coil voltages
  real*8              :: current_FB_fact  = 1.d0         !< Factor used for current feedback during the freeboundary equilibrium
  real*8, allocatable :: diag_coil_curr(:,:)

  type :: t_starwall_response
    integer :: file_version           = 9999
    integer :: n_bnd                  = -1
    integer :: nd_bez                 = -1
    integer :: ncoil                  = -1
    integer :: npot_w                 = -1
    integer :: n_w                    = -1
    integer :: ntri_w                 = -1
    integer :: n_tor                  = -1
    integer :: n_tor0                 = -1
    integer :: ntri_c                 = 0  !< Number of coil triangles
    integer :: n_pol_coils            = 0  !< Number of poloidal field coils
    integer :: n_rmp_coils            = 0  !< Number of RMP coils
    integer :: n_voltage_coils        = 0  !< Number of "voltage coils"
    integer :: n_diag_coils           = 0  !< Number of diagnostic coils
    integer :: ind_start_pol_coils    = -1 !< Index of first poloildal field coil
    integer :: ind_start_rmp_coils    = -1 !< Index of first RMP coil
    integer :: ind_start_voltage_coils= -1 !< Index of first voltage coil
    integer :: ind_start_diag_coils   = -1 !< Index of first diagnostic coil
    integer, allocatable :: jtri_c(:)      !< Number of triangles per coil
    real*8,  allocatable :: x_coil(:,:)      !< x-position of coil triangle nodes
    real*8,  allocatable :: y_coil(:,:)      !< y-position of coil triangle nodes
    real*8,  allocatable :: z_coil(:,:)      !< z-position of coil triangle nodes
    real*8,  allocatable :: phi_coil(:,:)    !< "Potential" at coil triangle nodes
    real*8,  allocatable :: eta_thin_coil(:) !< Thin wall resistivity of coil triangles
    real*8,  allocatable :: coil_resist(:)   !< Resistance of each coil
    real*8  :: eta_thin_w             = 1. !< Thin wall resistivity of wall triangles
    integer, allocatable :: i_tor(:)
    real*8,  allocatable :: d_yy(:)
    type(t_distrib_mat)  :: a_ye
    type(t_distrib_mat)  :: a_ey
    type(t_distrib_mat)  :: a_ee
    type(t_distrib_mat)  :: a_id
    type(t_distrib_mat)  :: a_nw
    type(t_distrib_mat)  :: s_ww
    type(t_distrib_mat)  :: s_ww_inv
    real*8,  allocatable :: xyzpot_w(:,:)
    integer, allocatable :: jpot_w(:,:)
  end type t_starwall_response
  
  type(t_starwall_response) :: sr              !< STARWALL response
  
  !> Coil current specification.
  type :: t_coil_curr_input
    real*8             :: current   = 0.d0  !< Current of the coil in Ampere*Turns
    real*8             :: pert      = 0.d0  !< Pert. of coil current in Ampere*Turns to speed-up VDE.
    real*8             :: pert_start_time  = 0.d0   !< Starting time of pert. of coil current in JOREK_time.
    real*8             :: pert_growth_time = 1.d-12 !< Ramp-up time of pert. of coil current in JOREK_time.
    character(len=256) :: curr_file = 'none'!< Ascii file with coil current time trace.
    real*8             :: time_shift    = 0.d0  !< Shift time of time trace.
    real*8             :: time_scale    = 1.d0  !< Scale time of time trace.
    real*8             :: curr_scale    = 1.d0  !< Scale amplitude of time trace.
    character(len=512) :: curr_expr = 'none'!< Analyt. expression for time trace (Python).
    real*8             :: max_time  = 0.d0  !< Evaluate analytical expression up to this time.
    integer            :: len       = 0     !< Evaluate analytical expression on this number of time points.
  end type t_coil_curr_input
  type :: t_coil_curr_time_trace
    integer            :: len      !< Number of points in numerical time trace.
    real*8, allocatable:: time(:)  !< time-values of numerical time trace
    real*8, allocatable:: curr(:)  !< current-values of numerical time trace
  end type t_coil_curr_time_trace
  integer, parameter              :: MAX_COILS = 299
  type(t_coil_curr_input), target :: diag_coils(MAX_COILS)
  type(t_coil_curr_input), target :: rmp_coils(MAX_COILS) 
  type(t_coil_curr_input), target :: voltage_coils(MAX_COILS)
  type(t_coil_curr_input), target :: pf_coils(MAX_COILS)
  type(t_coil_curr_time_trace)    :: coil_curr_time_trace(4*MAX_COILS)
  
  real*8 :: vert_FB_amp(MAX_COILS) = 0.d0 !< Allows to tune direction and magnitude of vertical feedback for each poloidal field coil.
  
  
  
  contains
  
  
  
  !> Set coil_curr_time_trace based on user input.
  subroutine set_coil_curr_time_trace()
    use profiles, only: readProf
    
    integer :: i, j, k, l, err
    character(len=60) :: s, filename
    real*8 :: r
    class(t_coil_curr_input), pointer :: coil_curr_input
    
    !--- Make sure that the main coil current vector (I_coils) is properly allocated
    if (sr%ncoil > 0) then
      n_coils = sr%ncoil
      if (.not. allocated(I_coils)) then 
        allocate(I_coils(n_coils))
        I_coils = 0.d0
      endif
    endif
        
    do i = 1, sr%ncoil
            
      ! --- Which coil are we processing?
      if ( (i>=sr%ind_start_pol_coils) .and. (i<sr%ind_start_pol_coils+sr%n_pol_coils) ) then
        j = i - sr%ind_start_pol_coils + 1
        coil_curr_input => pf_coils(j)
      else if ( (i>=sr%ind_start_rmp_coils) .and. (i<sr%ind_start_rmp_coils+sr%n_rmp_coils) ) then
        j = i - sr%ind_start_rmp_coils + 1
        coil_curr_input => rmp_coils(j)
      else if ( (i>=sr%ind_start_diag_coils) .and. (i<sr%ind_start_diag_coils+sr%n_diag_coils) ) then
        j = i - sr%ind_start_diag_coils + 1
        coil_curr_input => diag_coils(j)
      else if ( (i>=sr%ind_start_voltage_coils) .and. (i<sr%ind_start_voltage_coils+sr%n_voltage_coils) ) then
        j = i - sr%ind_start_voltage_coils + 1
        coil_curr_input => voltage_coils(j)
      end if
      
      ! --- Set up numerical coil current time trace based on the input type...
      if ( coil_curr_input%curr_file /= 'none' ) then ! ... numerical input by ascii file
        
        call readProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, &
          coil_curr_time_trace(i)%len, coil_curr_input%curr_file)
        
      else if ( coil_curr_input%curr_expr /= 'none' ) then ! ... analytical Python expression
        
        ! --- Python script
        call random_seed()
        err = 1
        do while ( err /= 0 )
          call random_number(r)
          l = r * 99999999
          write(s,*) l
          filename='./jorek_curr_expr_'//trim(adjustl(s))//'.py'
          open(42, file=trim(filename), status='new', iostat=err)
        end do
        111 format(2a)
        112 format(a,i16)
        113 format(a,es25.16)
        write(42,111) 'from math import *'
        write(42,111) 'def f(t):'
        write(42,111) '  return ', trim(coil_curr_input%curr_expr)
        write(42,112) 'len=', coil_curr_input%len
        write(42,113) 'tmin=', - coil_curr_input%time_shift / coil_curr_input%time_scale - 1.d-12
        write(42,113) 'tmax=', ( coil_curr_input%max_time - coil_curr_input%time_shift ) / coil_curr_input%time_scale + 1.d-12
        write(42,111) 'for x in range(1,len):'
        write(42,111) '  t=tmin+(x-1)/float(len-1)*(tmax-tmin)'
        write(42,111) '  s = "%25.16e"%t'
        write(42,111) '  s += "%25.16e"%f(t)'
        write(42,111) '  print(s)'
        close(42)
        
        ! --- Call Python
        call system('python ./jorek_curr_expr_'//trim(adjustl(s))//'.py > ./jorek_curr_expr_'//trim(adjustl(s))//'.dat')
        
        ! --- Read the result
        call readProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, &
          coil_curr_time_trace(i)%len, './jorek_curr_expr_'//trim(adjustl(s))//'.dat')
        
        ! --- Delete temporary files
        call system('rm ./jorek_curr_expr_'//trim(adjustl(s))//'.py ./jorek_curr_expr_'//trim(adjustl(s))//'.dat')
        
      else ! ... just a value plus a perturbation
        
        if (allocated(coil_curr_time_trace(i)%time)) deallocate(coil_curr_time_trace(i)%time)
        if (allocated(coil_curr_time_trace(i)%curr)) deallocate(coil_curr_time_trace(i)%curr)
        
        allocate(coil_curr_time_trace(i)%time(4) )
        allocate(coil_curr_time_trace(i)%curr(4) )
        coil_curr_time_trace(i)%len = 4
        
        coil_curr_time_trace(i)%time(1)   = -1.d12
        coil_curr_time_trace(i)%time(2)   = coil_curr_input%pert_start_time
        coil_curr_time_trace(i)%time(3)   = coil_curr_input%pert_start_time + coil_curr_input%pert_growth_time
        coil_curr_time_trace(i)%time(4)   = 1.d12
        
        coil_curr_time_trace(i)%curr(1:2) = coil_curr_input%current
        coil_curr_time_trace(i)%curr(3:4) = coil_curr_input%current + coil_curr_input%pert
        
      end if
      
      coil_curr_time_trace(i)%time = (coil_curr_time_trace(i)%time * coil_curr_input%time_scale) &
        + coil_curr_input%time_shift
      coil_curr_time_trace(i)%curr = coil_curr_time_trace(i)%curr * coil_curr_input%curr_scale
      
      if ( coil_curr_time_trace(i)%time(1) > 0.d0 ) then
        write(*,*) 'ERROR: A coil current time trace does not start at time 0. Maybe you have time_shift wrong?', i
        stop
      end if
      
    end do
    
  end subroutine set_coil_curr_time_trace
  
  
  
  !> Check whether the user has supplied the right coil current time trace input.
  subroutine check_coil_curr_time_trace_input(coils_number)
  
    integer, intent(in) :: coils_number
    
    integer :: i, j
    logical :: changed_by_user
    class(t_coil_curr_input), pointer :: coil_curr_input
    
    do i = 1, 4
      do j = 1, MAX_COILS
        
        if ( i == 1 ) then
          coil_curr_input => pf_coils(j)
        else if ( i == 2 ) then
          coil_curr_input => rmp_coils(j)
        else if ( i == 3 ) then
          coil_curr_input => diag_coils(j)
        else if ( i == 4 ) then
          coil_curr_input => voltage_coils(j)
        end if
        
        changed_by_user = ( (coil_curr_input%current    /= 0.d0  ) .or. (coil_curr_input%pert       /= 0.d0)  &
                       .or. (coil_curr_input%curr_file  /= 'none') .or. (coil_curr_input%time_shift /= 0.d0)  &
                       .or. (coil_curr_input%time_scale /= 1.d0  ) .or. (coil_curr_input%curr_scale /= 1.d0)  &
                       .or. (coil_curr_input%curr_expr  /= 'none') .or. (coil_curr_input%max_time   /= 0.d0)  &
                       .or. (coil_curr_input%len        /= 0     ) .or. (vert_FB_amp(j)             /= 0.d0)  )

        if ( (j > coils_number) .and. changed_by_user ) then
          write(*,*) 'ERROR: Coil current input has been provided for a coil not existing.', i, j
          stop
        end if
        
        if ( (coil_curr_input%curr_expr /= 'none') .and. ( (coil_curr_input%len <= 0) .or. (coil_curr_input%max_time <= 0.d0) ) ) then
          write(*,*) 'ERROR: When prescribing a coil current via %curr_expr, you need to specify also %len and %max_time.'
          stop
        end if
        
      end do
    end do
  end subroutine check_coil_curr_time_trace_input
  
  
  
  !> Preset freeboundary related input parameters to reasonable default values.
  subroutine vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(out)   :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Preset namelist input parameters.
    freeboundary_equil   = .false.
    starwall_equil_coils = .false.
    freeboundary         = .false.
    resistive_wall       = .true.
    wall_resistivity     = 0.d0
    wall_resistivity_fact= 1.d0
        
    current_ref          = 1.d22
    FB_Ip_position       = 0.2d0
    FB_Ip_integral       = 0.01d0
    n_feedback_current   = 2
        
    Z_axis_ref           = 1.d22
    FB_Zaxis_position    = 1.d0
    FB_Zaxis_derivative  = 0.d0
    FB_Zaxis_integral    = 0.d0
    n_feedback_vertical  = 1
    start_VFB            = 10
    
    n_iter_freeb         = 900
    
    PF_pert_start_time   = 1.d99
    psi_offset_freeb     = 0.d0
    
  end subroutine vacuum_preset
  
  
  
  !> Initialize vacuum, ensure consistency of input parameters.
  subroutine vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)
    
    integer, intent(in)    :: my_id
    logical, intent(inout) :: freeboundary_equil, freeboundary, resistive_wall
    
    ! --- Make input parameters consistent.
    freeboundary   = freeboundary .or. freeboundary_equil
    resistive_wall = freeboundary .and. resistive_wall
        
    ! --- Initialize some variables.
    sr%n_bnd  = 0
    sr%nd_bez = 0
    sr%ncoil  = 0
    sr%npot_w = 0
    sr%n_w    = 0
    sr%ntri_w = 0
    sr%n_tor  = 0
    sr%n_tor0 = 0
    
    ! --- Switch on terms on the RHS of current equation definition when using free-boundary
    freeb_fact = 0.d0
    if ( freeboundary ) freeb_fact = 1.d0
    
    if ( (my_id == 0) .and. (sum(pf_coils%pert) > 0) .and. (PF_pert_start_time>1.d30) ) then
       write(*,*) 'WARNING: Poloidal field coil perturbation pf_coils%pert has been set by the user, but will not be applied since PF_pert_start_time was not set to a reasonable value.'
    end if
    
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
    do i = 1, sr%n_tor
      is_freebound = is_freebound .or. ( i_tor == sr%i_tor(i) )
    end do
    
    ! --- Free boundary conditions only for certain variables
    is_freebound = is_freebound .and. ( (i_var == ivar_j) .or. (i_var == ivar_psi) )
    
  end function is_freebound
  
  
  
  !> Import some vacuum-related data from the restart file
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_restart_vacuum(file_handle, freeboundary, resistive_wall)
    
    use phys_module, only: t_start
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    logical :: resistive_wall_rst
    integer :: ierr
    
    if ( freeboundary ) then
      
      read(file_handle, iostat=ierr) resistive_wall_rst
      if ( ierr /= 0 ) then
        write(*,*) 'WARNING: Restarting a simulation with freeboundary=.t. which was run with'
        write(*,*) '  freeboundary=.f. so far.'
        return
      else if ( resistive_wall .neqv. resistive_wall_rst ) then
        write(*,*) 'ERROR: It is currently not possible to restart a JOREK simulation with a'
        write(*,*) '  modified setting for resistive_wall.'
        stop
      end if
      
      if ( resistive_wall ) then
        
        read(file_handle) n_wall_curr, n_dof_starwall
        
        if ( allocated(wall_curr) ) deallocate(wall_curr)
        allocate( wall_curr(n_wall_curr) )
        read(file_handle) wall_curr(:)
        
        if ( allocated(dwall_curr) ) deallocate(dwall_curr)
        allocate( dwall_curr(n_wall_curr) )
        read(file_handle) dwall_curr(:)
        
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
        old_dpsibnd_vec = 0.d0 !###
        
        if ( vacuum_debug .and. resistive_wall ) then
          write(*,*) 'DEBUG: Checksums'
          write(*,*) 'wall_curr', sum(abs(wall_curr))
          write(*,*) 'dwall_curr', sum(abs(dwall_curr))
          write(*,*) 'END: Checksums'
        end if
        
        wall_curr_initialized = .true.
        
      end if
      
      read(file_handle) current_FB_fact
      read(file_handle) n_coils
      if ( n_coils /= 0 ) then
        if ( allocated(I_coils) ) deallocate(I_coils)
        allocate( I_coils(n_coils) )
        read(file_handle) I_coils(:)
      end if
      
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine import_restart_vacuum


  !> Import some vacuum-related data from the HDF5 restart file
  !!
  !! @todo Does not work currently if variable freeboundary is changed between export and import!
  subroutine import_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
    
    use phys_module, only: t_start, nstep, index_start
#ifdef USE_HDF5
    use hdf5
    use hdf5_io_module
    
    ! --- Routine parameters
    integer(HID_T), intent(in) :: file_id
#else
    integer,        intent(in) :: file_id
#endif
    logical,        intent(in) :: freeboundary
    logical,        intent(in) :: resistive_wall
    
#ifdef USE_HDF5
    ! --- Local variables
    logical   :: freeboundary_rst, resistive_wall_rst
    integer   :: n_diag_coil = 0
    character :: t_freeboundary, t_resistive_wall
    real*8, allocatable :: t_diag_coil_curr(:,:)
    
    call HDF5_char_reading(file_id,t_freeboundary,"freeboundary")
    freeboundary_rst = (t_freeboundary == "T")
    
    if ( freeboundary ) then
      
      if ( .not. freeboundary_rst ) then
        write(*,*) 'WARNING: Restarting a simulation with freeboundary=.t. which was run with'
        write(*,*) '  freeboundary=.f. so far.'
      end if
      
      if ( freeboundary_rst ) then
        call HDF5_char_reading(file_id,t_resistive_wall,"resistive_wall")
        resistive_wall_rst = (t_resistive_wall == "T")
      else
        resistive_wall_rst = .false.
      end if
      
      if ( resistive_wall .and. resistive_wall_rst ) then
        
        call HDF5_integer_reading(file_id,n_wall_curr,"n_wall_curr")
        call HDF5_integer_reading(file_id,n_dof_starwall,"n_dof_starwall")
        
        if ( allocated(wall_curr) ) deallocate(wall_curr)
        allocate( wall_curr(n_wall_curr) )
        call HDF5_array1D_reading(file_id,wall_curr,"wall_curr")
        
        if ( allocated(dwall_curr) ) deallocate(dwall_curr)
        allocate( dwall_curr(n_wall_curr) )
        call HDF5_array1D_reading(file_id,dwall_curr,"dwall_curr")
        
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
        old_dpsibnd_vec(:) = 0.d0
        call HDF5_array1D_reading(file_id,old_dpsibnd_vec,"old_dpsibnd_vec")
        
        if ( index_start > 1 ) then
          if ( allocated(diag_coil_curr) ) deallocate(diag_coil_curr)
          call HDF5_integer_reading(file_id,n_diag_coil,"n_diag_coil")
          allocate( t_diag_coil_curr(index_start,n_diag_coil) )
          call HDF5_array2D_reading(file_id,t_diag_coil_curr,"diag_coil_curr")
          allocate( diag_coil_curr(index_start+nstep,n_diag_coil) )
          diag_coil_curr(1:index_start,:) = t_diag_coil_curr(1:index_start,:)
          deallocate(t_diag_coil_curr)
        end if
        
        if ( vacuum_debug .and. resistive_wall ) then
          write(*,*) 'DEBUG: Checksums'
          write(*,*) 'wall_curr', sum(abs(wall_curr))
          write(*,*) 'dwall_curr', sum(abs(dwall_curr))
          write(*,*) 'END: Checksums'
        end if
        
        wall_curr_initialized = .true.
        
      else if ( resistive_wall ) then
        
        write(*,*) 'WARNING: Continuing a JOREK simulation with resistive_wall=.f.'
        write(*,*) '   now with resistive_wall=.t. Wall currents will be initialized to zero.'
        wall_curr_initialized = .false.
        
      end if
      
      call HDF5_real_reading(file_id,current_FB_fact,'current_FB_fact')
      call HDF5_integer_reading(file_id,n_coils,"n_coils")
      if ( n_coils /= 0 ) then
        if ( allocated(I_coils) ) deallocate(I_coils)
        allocate( I_coils(n_coils) )
        call HDF5_array1D_reading(file_id,I_coils,"I_coils")
      end if
      
    end if
     
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if

#endif

  end subroutine import_HDF5_restart_vacuum
  
  
  !> Export some vacuum-related data to the restart file
  subroutine export_restart_vacuum(file_handle, freeboundary, resistive_wall)
    
    ! --- Routine parameters
    integer, intent(in) :: file_handle
    logical, intent(in) :: freeboundary
    logical, intent(in) :: resistive_wall
    
    if ( freeboundary ) then
      
      write(file_handle) resistive_wall
      if ( resistive_wall ) then
        
        if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or.                    &
          (.not. allocated(old_dpsibnd_vec)) ) then
          write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
          stop
        end if
        
        write(file_handle) n_wall_curr, n_dof_starwall
        
        write(file_handle) wall_curr(:)
        
        write(file_handle) dwall_curr(:)
        
      end if
      
      write(file_handle) current_FB_fact
      
      if ( (n_coils/=0) .and. (.not. allocated(I_coils)) ) then
        write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: I_coils not allocated.'
        stop
      end if
    
      write(file_handle) n_coils
      if ( n_coils /= 0 ) write(file_handle) I_coils(:)
      
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
    
  end subroutine export_restart_vacuum
  
  
  !> Export some vacuum-related data to the HDF5 restart file
  subroutine export_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
    
    use phys_module, only: index_now
#ifdef USE_HDF5
    use hdf5
    use hdf5_io_module

    ! --- Routine parameters
    integer(HID_T), intent(in) :: file_id
#else
    integer,        intent(in) :: file_id
#endif
    logical,        intent(in) :: freeboundary
    logical,        intent(in) :: resistive_wall
    
#ifdef USE_HDF5
    ! --- Local variables
    character           :: t_freeboundary, t_resistive_wall
    real*8, allocatable :: t_diag_coil_curr(:,:)

    t_freeboundary = "F"
    if (freeboundary) t_freeboundary = "T"
    call HDF5_char_saving(file_id,t_freeboundary,"freeboundary"//char(0))
    
    if ( freeboundary ) then
    
      t_resistive_wall = "F"
      if (resistive_wall) t_resistive_wall = "T"
      call HDF5_char_saving(file_id,t_resistive_wall,"resistive_wall"//char(0))

      if ( resistive_wall ) then
        if ( (.not. allocated(wall_curr)) .or. (.not. allocated(dwall_curr)) .or. &
           (.not. allocated(old_dpsibnd_vec)) ) then
          write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: Arrays not allocated.'
          stop
        end if
        call HDF5_integer_saving(file_id,n_wall_curr,"n_wall_curr"//char(0))
        call HDF5_integer_saving(file_id,n_dof_starwall,"n_dof_starwall"//char(0))
        call HDF5_array1D_saving(file_id,wall_curr,n_wall_curr,"wall_curr"//char(0))
        call HDF5_array1D_saving(file_id,dwall_curr,n_wall_curr,"dwall_curr"//char(0))
        call HDF5_integer_saving(file_id,sr%n_diag_coils,"n_diag_coil"//char(0))
        
        if ( (index_now > 0) .and. (sr%n_diag_coils > 0) ) then
          allocate(t_diag_coil_curr(index_now,sr%n_diag_coils))
          t_diag_coil_curr(1:index_now,:) = diag_coil_curr(1:index_now,:)
          call HDF5_array2D_saving(file_id,t_diag_coil_curr,index_now,sr%n_diag_coils,"diag_coil_curr"//char(0))
          deallocate(t_diag_coil_curr)
        end if
      end if

      call HDF5_array1D_saving(file_id,old_dpsibnd_vec,n_dof_starwall,"old_dpsibnd_vec"//char(0))
      
      call HDF5_real_saving(file_id,current_FB_fact,'current_FB_fact'//char(0))
      if ( (n_coils/=0) .and. (.not. allocated(I_coils)) )  then
        write(*,*) 'ERROR in mod_vacuum.f90:export_restart_vacuum: I_coils not allocated.'
        stop
      end if
      call HDF5_integer_saving(file_id,n_coils,"n_coils"//char(0))
      if ( n_coils /= 0 ) call HDF5_array1D_saving(file_id,I_coils,n_coils,"I_coils"//char(0))
    end if
    
    if ( vacuum_debug .and. resistive_wall ) then
      write(*,*) 'DEBUG: Checksums'
      if ( allocated(wall_curr)  ) write(*,*) 'wall_curr ', sum(abs(wall_curr))
      if ( allocated(dwall_curr) ) write(*,*) 'dwall_curr', sum(abs(dwall_curr))
      write(*,*) 'END: Checksums'
    end if
#endif
  
  end subroutine export_HDF5_restart_vacuum
  
  
  !> Broadcast vacuum information between MPI processes
  subroutine broadcast_vacuum(my_id, resistive_wall)
    use mpi_mod
    implicit none
    
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    integer :: ierr, sz(2)
    
    call MPI_BCAST(n_dof_starwall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    if ( resistive_wall ) then
      
      call MPI_BCAST(n_wall_curr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      if ( n_wall_curr == 0 ) return
      
      if ( my_id == 0 ) sz(:) = (/ size(diag_coil_curr,1), size(diag_coil_curr,2) /)
      call MPI_BCAST(sz,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      if ( my_id /= 0 ) then
        if ( allocated(wall_curr) ) deallocate(wall_curr)
        allocate( wall_curr(n_wall_curr) )
        
        if ( allocated(dwall_curr) ) deallocate(dwall_curr)
        allocate( dwall_curr(n_wall_curr) )
        
        if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
        allocate( old_dpsibnd_vec(n_dof_starwall) )
        
        if ( allocated(diag_coil_curr) ) deallocate(diag_coil_curr)
        if ( minval(sz) > 0 ) allocate( diag_coil_curr(sz(1),sz(2)) )
      else
      end if
      
      call MPI_BCAST(wall_curr,n_wall_curr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dwall_curr,n_wall_curr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(old_dpsibnd_vec,n_dof_starwall,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
      call MPI_BCAST(wall_curr_initialized,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      if ( minval(sz) > 0 ) call MPI_BCAST(diag_coil_curr,sz(1)*sz(2),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      
    end if
    
    call MPI_BCAST(current_FB_fact,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(freeb_fact,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
  end subroutine broadcast_vacuum
  
  
  
end module vacuum

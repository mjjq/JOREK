!> Module for execution of user commands
!! (used by jorek2_postproc)
!!
!! @todo Correct treatment of cases without xpoint or double xpoints
module exec_commands
  
  use parameters,        only: n_var, n_tor, n_order, n_vertex_max, variable_names, jorek_model
  use phys_module,       only: mode, t_start, xpoint, xcase, index_start, F0, ZK_perp, ZK_par,     &
    LOWER_XPOINT, UPPER_XPOINT, DOUBLE_NULL
  use data_structure,    only: type_node_list, type_element_list, type_surface_list, type_surface, &
    type_node, type_element
  use nodes_elements,    only: node_list, element_list
  use parse_commands,    only: type_command, print_command
  use settings,          only: set_setting, get_setting, get_int_setting, print_settings
  use convert_character, only: to_float, to_int, get_variable_number, lower_case
  use postproc_help,     only: general_help, specific_help
  use basis_at_gaussian, only: H, H_s, H_t, wgauss, n_gauss
  use domains,           only: in_private
  
  
  
  implicit none
  
  
  
  !> Routine interface to check if the number of command parameters is
  !! * equal to the given value (check_param_count1)
  !! * or equal to one of the given values (check_param_count2)
  interface check_param_count
    module procedure check_param_count1
    module procedure check_param_count2
  end interface check_param_count
  
  integer, parameter :: NORMAL_MODE = 1 !< Normal mode
  integer, parameter :: LOOP_MODE   = 2 !< Mode started by 'for' ended by 'done' commands
  integer :: exec_mode = NORMAL_MODE    !< Current operation mode (NORMAL_MODE or LOOP_MODE)
  integer :: loop_min_step              !< Smallest timestep index for current loop
  integer :: loop_max_step              !< Largest timestep index for current loop
  
  integer, parameter :: MAX_QUEUE_LENGTH  = 16384        !< Maximum length of command queue
  integer            :: n_queued_commands = 0            !< Number of commands in the queue
  type(type_command) :: command_queue(MAX_QUEUE_LENGTH)  !< Queued commands
  
  logical :: input_loaded  = .false.    !< Has an input file already been loaded?
  logical :: step_imported = .false.    !< Has a restart file already been imported?
  
  !> @name Positions of x-point and magnetic axis (determined by load_step)
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  real*8  :: psi_bnd, psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
  integer :: i_elm_axis, i_elm_xpoint(2)
  
  
  
  private
  public exec_command, general_help, specific_help

  
  
  save
  
  
  
  contains
  
  
  
  !> Checks if the number of parameters is correct (first version); use interface check_param_count
  subroutine check_param_count1(command, n_params, error)
    ! --- Routine parameters
    type(type_command), intent(in)  :: command  !< Command to check
    integer,            intent(in)  :: n_params !<  Required number of parameters to the command
    integer,            intent(out) :: error    !< Error flag
    
    call check_param_count2(command, (/n_params/), error)
    
  end subroutine check_param_count1
  
  
  
  !> Checks if the number of parameters is correct (second version); use interface check_param_count
  subroutine check_param_count2(command, n_params, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command      !< Command to check
    integer,            intent(in)  :: n_params(:)  !< Required number of parameters to the command
    integer,            intent(out) :: error        !< Error flag
    
    ! --- Local variables
    integer :: i
    
    error = 0
    
    do i = 1, size(n_params)
      if ( command%n_options - 1 == n_params(i) ) return
    end do
    write(*,'(1x,a,a,a,99i3)') 'ERROR (command ', trim(command%option(1)), &
      '): Required number of parameters is one of:', n_params(:)
    
    call specific_help(command%option(1))
    error = 1
    
  end subroutine check_param_count2
  
  
  
  !> Determines the value of variable ivar at position (R,Z,phi)
  real*8 function variable_value(ivar, R, Z, phi, error)
    
    ! --- Routine parameters
    integer, intent(in) :: ivar      !< Number of the variable
    real*8,  intent(in) :: R         !< R-position to evaluate the variable at
    real*8,  intent(in) :: Z         !< Z-position to evaluate the variable at
    real*8,  intent(in) :: phi       !< phi-position to evaluate the variable at   
    integer, intent(out):: error     !< Error flag
    
    ! --- Local variables
    integer :: ielm, i_harm
    real*8  :: R_out, Z_out, s, t, P, P_s, P_t, P_st, P_ss, P_tt
    
    error = 0
    variable_value = 0.d0
    do i_harm = 1, n_tor
    
      call find_RZ(node_list, element_list, R, Z, R_out, Z_out, ielm, s, t, error)
      if ( error /= 0 ) return
      call interp (node_list, element_list, ielm, ivar, i_harm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
      
      if ( i_harm == 1 ) then
         variable_value = variable_value + P
      else if ( mod(i_harm,2) == 0 ) then
         variable_value = variable_value + P * cos( mode(i_harm) * phi )
      else
         variable_value = variable_value + P * sin( mode(i_harm) * phi )
      end if
      
    end do
    
  end function variable_value
  
  
  
  !> Determines the value of toroidal mode i_harm of variable ivar at position (R,Z)
  real*8 function variable_mode_value(ivar, R, Z, i_harm, error)
   
    ! --- Routine parameters
    integer, intent(in) :: ivar      !< Number of the variable
    real*8,  intent(in) :: R         !< R-position to evaluate the variable at
    real*8,  intent(in) :: Z         !< Z-position to evaluate the variable at
    integer, intent(in) :: i_harm
    integer, intent(out):: error     !< Error flag
    
    ! --- Local variables
    integer :: ielm
    real*8  :: R_out, Z_out, s, t, P_s, P_t, P_st, P_ss, P_tt
   
    error = 0
    
    call find_RZ(node_list, element_list, R, Z, R_out, Z_out, ielm, s, t, error)
    if ( error /= 0 ) return
    call interp (node_list, element_list, ielm, ivar, i_harm, s, t, variable_mode_value, P_s, P_t, &
      P_st, P_ss, P_tt)
    
  end function variable_mode_value
  
  
  
  !> Execute a command
  recursive subroutine exec_command(command, first_step, file_handle, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    error = 0
    
    if ( get_setting('debug',error) == 'true' ) then
      write(*,'(a)') 'Exec_command was called with:'
      call print_command(command)
      write(*,*)
    end if
    
    ! --- In normal mode, some commands are directly executed
    if ( exec_mode == NORMAL_MODE ) then
      
      if ( ( .not. input_loaded ) .and. ( trim(command%option(1)) == 'for' ) ) then
        write(*,*) 'ERROR: No namelist input file loaded.'
        call specific_help('namelist')
        error = 1
        return
      end if
      
      select case ( trim(command%option(1)) )
        case ( 'average' )
          call average(command, first_step, file_handle, error)
        case ( 'axis' )
          call magnetic_axis(command, first_step, file_handle, error)
        case ( 'global_parameters' )
          call global_parameters(command, first_step, file_handle, error)
        case ( 'fluxsurfaces' )
          call fluxsurfaces(command, file_handle, error)
        case ( 'for' )
          call loop_start(command, error)
        case ( 'heatfluxpattern' )
          call heatfluxpattern(command, file_handle, error)
        case ( 'help' )
          call help(command, error)
        case ( 'line' )
          call line(command, first_step, file_handle, error)
        case ( 'namelist' )
          call load_namelist(command, file_handle, error)
        case ( 'params' )
          call log_parameters(0)
        case ( 'point' )
          call point(command, first_step, file_handle, error)
        case ( 'set' )
          call set(command, error)
        case ( 'settings' )
          call print_settings()
        case ( 'timesteps' )
          call timesteps 
        case ( 'volume' )
          call plasma_volume(command, first_step, file_handle, error)
        case default
          write(*,*) 'Command "', trim(command%option(1)), '" does not exist'
          call general_help() 
      end select
      
    ! --- In loop mode, commands are queued and afterwards executed for each timestep separately
    else if ( exec_mode == LOOP_MODE ) then
      
      select case ( trim(command%option(1)) )
        case ( 'axis', 'average', 'fluxsurfaces', 'heatfluxpattern', 'line', 'point', 'volume',    &
          'global_parameters' )
          call add_to_command_queue(command, error)
        case ( 'help' )
          call help(command, error)
        case ( 'done' )
          call loop_end(command, error)
        case default
          write(*,*) 'Command "', trim(command%option(1)), '" unknown/invalid inside a loop'
          call general_help() 
      end select
        
    end if
    
  end subroutine exec_command
  
  
  
  !> Loads a time step from a restart file if the restart file exists
  subroutine load_step(istep, error)
    
    ! --- Routine parameters
    integer,            intent(in)  :: istep !< Load this time step
    integer,            intent(out) :: error !< Error flag
    
    character(len=64) :: file_name
    logical           :: file_exists
    integer           :: ifail
    
    error = 0
    
    write(file_name,'(a,i5.5,a)') 'jorek', istep, '.rst'
    
    inquire(file=file_name, exist=file_exists)
    if ( .not. file_exists ) then
      error = 1
      return
    end if
    
    write(*,*)
    write(*,'(a,i5.5,a)') '#################### TIME STEP ', istep, ' ####################'
    write(*,*)
    
    ! --- Load the restart file
    call import_restart(node_list, element_list, trim(file_name), error)
    if ( error /= 0 ) return
    step_imported = .true.
    
    ! --- Locate magnetic axis and X-point.
    call find_axis(0, node_list, element_list, psi_axis, R_axis, Z_axis, i_elm_axis, s_axis,       &
      t_axis, ifail)
    if ( xpoint ) then
      call find_xpoint(0, node_list, element_list, psi_xpoint, R_xpoint, Z_xpoint, i_elm_xpoint,     &
        s_xpoint, t_xpoint, xcase, ifail)
      if ( (xcase == DOUBLE_NULL) ) then
        psi_bnd = minval( psi_xpoint(:) )
      else
        psi_bnd = psi_xpoint(1)
      end if
    else
      psi_bnd = 0.d0
    end if
    ! ### handle cases without x-point or with double x-point ###
    
  end subroutine load_step
  
  
  
  !> Is called at the start of a for loop.
  !!
  !! It switches the module behaviour (exec_mode) to LOOP_MODE, i.e., all commands inside the for
  !! loop are not executed immediately but collected in a command queue. All commands of the
  !! command queue are executed for each time step after the 'done' command of the for loop
  !! has been entered.
  subroutine loop_start(command, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: error   !< Error flag
    
    error = 0
    call check_param_count(command, (/3,5/), error)
    if ( error /= 0 ) return
    
    if ( trim(command%option(2)) /= 'step' ) then
      write(*,*) 'ERROR: Wrong syntax for "for" statement.'
      call specific_help('for')
      return
    end if
    
    loop_min_step = to_int(command%option(3),error)
    if ( error /= 0 ) then
      write(*,*) 'ERROR: Wrong syntax for "for" statement.'
      call specific_help('for')
      return
    end if
    
    if ( command%n_options == 4 ) then ! for step <XXX> do
      if ( trim(command%option(4)) /= 'do' ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
      loop_max_step = loop_min_step
    else ! for step <XXX> to <YYY> do
      if ( trim(command%option(4)) /= 'to' ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
      loop_max_step = to_int(command%option(5),error)
      if ( error /= 0 ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
      if ( trim(command%option(6)) /= 'do' ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
    end if
    
    exec_mode = LOOP_MODE
    
  end subroutine loop_start
  
  
  
  !> Is called upon the 'done' command of a for loop and executes all commands enclosed in the loop.
  subroutine loop_end(command, error)
    
    include 'omp_lib.h'
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: error   !< Error flag
    
    ! --- Local variables
    integer :: jcmd, istep, load_error
    integer :: thread_id, file_handle
    logical :: first_step ! Is true for the first timestep loaded in the for-loop
    
    error = 0
    exec_mode = NORMAL_MODE    
    
    if ( n_queued_commands == 0 ) then
      write(*,*) 'WARNING: No commands in for-loop.'
      return
    end if
    
    first_step = .true.
    do istep = loop_min_step, loop_max_step
      call load_step(istep, load_error)
      if ( load_error /= 0 ) cycle
      
!$omp parallel do default(shared) private(error,thread_id,file_handle)
      do jcmd = 1, n_queued_commands
        
        thread_id   = omp_get_thread_num()
        file_handle = 17 + thread_id
        
        call exec_command(command_queue(jcmd), first_step, file_handle, error)  
        if ( error /= 0 ) then
          write(*,*) 'ERROR executing the following command (ignoring it):'
          call print_command(command_queue(jcmd))
          error = 0
        end if
        
      end do
!$omp end parallel do
      
      first_step = .false.
    end do
    
    if ( first_step ) then
      write(*,'(a,i5.5,a,i5.5,a)') 'WARNING: There were no restart files for steps ',              &
        loop_min_step, ' to ', loop_max_step, '.'
    end if
    
    n_queued_commands = 0
    
  end subroutine loop_end
  
  
  
  !> Add a command to the command queue (for loop)
  subroutine add_to_command_queue(command, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: error   !< Error flag
    
    error = 0
    
    if ( n_queued_commands + 1 > MAX_QUEUE_LENGTH) then
      write(*,*) 'ERROR: Too many commands in the command queue of the for loop.'
      error = 1
      return
    end if
    
    n_queued_commands = n_queued_commands + 1
    command_queue(n_queued_commands) = command
    
  end subroutine add_to_command_queue
  
  
  
  !> Retrieve a command from the command queue (currently not used!)
  subroutine get_from_command_queue(command, error)
    
    ! --- Routine parameters
    type(type_command), intent(out)    :: command !< Command to be executed
    integer,            intent(out)    :: error   !< Error flag
    
    error = 0
    
    if ( n_queued_commands < 1 ) then
      write(*,*) 'ERROR: Cannot get a command from an empty queue.'
      error = 1
      return
    end if
    
    n_queued_commands = n_queued_commands - 1
    command = command_queue(1)
    command_queue(1:n_queued_commands) = command_queue(2:n_queued_commands+1)
    
  end subroutine get_from_command_queue
  
  
  
  !> Implements the 'fluxsurfaces' command: Writes out the flux surfaces
  recursive subroutine fluxsurfaces(command, file_handle, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command     !< Command to be executed
    integer,            intent(in)     :: file_handle !< File handle
    integer,            intent(out)    :: error       !< Error flag
    
    ! --- Local variables
    integer                  :: i_elm
    integer                  :: i, j, k, nplot, ip
    real*8                   :: ss1, dss1, ss2, dss2, tt1, dtt1, tt2, dtt2
    real*8                   :: u, si, dsi, ti, dti
    real*8                   :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    character(len=1024)      :: filename
    type (type_surface_list) :: surface_list
    
    error = 0
    
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ! --- Open file
    write(filename,'(a,i5.5,a)') 'fluxsurfaces_step', index_start, '.dat'
    open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
    if (error /= 0) then
      write(*,*) 'ERROR in routine fluxsurfaces: Creating file "', trim(filename), '" failed.'
      return
    end if
    
    nplot = 11 ! Number of points in each flux surface segment
    surface_list%n_psi = get_int_setting('surfaces', error) ! Number of flux surfaces
    allocate (surface_list%psi_values(surface_list%n_psi))
    
    !### TODO: Find psi_max and use it ###
    !### TODO: use find_flux_surfaces2 instead ###
    
    do k = 1, surface_list%n_psi
      surface_list%psi_values(k) = psi_axis + ( psi_bnd - psi_axis) * 1.2d0 * real(k-1)            &
        / real(surface_list%n_psi-1)
    end do
    
    ! --- Find flux surfaces
    call find_flux_surfaces(xpoint,xcase,node_list,element_list,surface_list)
    
    ! --- Loop over all flux surfaces
    do i=1, surface_list%n_psi
      
      ! --- Loop over all segments of this flux surface
      do j=1,surface_list%flux_surfaces(i)%n_pieces
        
        ! --- Bezier element, in which the current flux surface segment is located
        i_elm = surface_list%flux_surfaces(i)%elm(j)
        
        ss1  = surface_list%flux_surfaces(i)%s(1,j)
        dss1 = surface_list%flux_surfaces(i)%s(2,j)
        ss2  = surface_list%flux_surfaces(i)%s(3,j)
        dss2 = surface_list%flux_surfaces(i)%s(4,j)
        
        tt1  = surface_list%flux_surfaces(i)%t(1,j)
        dtt1 = surface_list%flux_surfaces(i)%t(2,j)
        tt2  = surface_list%flux_surfaces(i)%t(3,j)
        dtt2 = surface_list%flux_surfaces(i)%t(4,j)
        
        ! --- Loop over nplot points in a flux surface segment
        do ip = 1, nplot
          u = -1. + 2.*float(ip-1)/float(nplot-1)
          
          ! --- Determine s and t values of the current point inside element i_elm
          call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
          call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
          
          ! --- Determine (R,Z)-coordinates of the current point on the current flux surface
          call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt, &
            Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
            
          ! --- Write out the (R,Z)-coordinates
          write(file_handle,'(2ES16.7)') R, Z
        end do
        
        write(file_handle,*)
        write(file_handle,*)

      end do
      
    end do

    close (unit=file_handle)
    deallocate (surface_list%psi_values)
    
  end subroutine fluxsurfaces
  
  
  
  !> Returns (R,Z)-coordinates of points on the n_psi flux surfaces specified by the psi_values
  recursive subroutine find_flux_surfaces2(n_psi, psi_values, n_pts, n_seg, R, Z, error)
    
    ! --- Routine parameters
    integer,                  intent(in)        :: n_psi              !< Number of flux surfaces
    real*8,                   intent(in)        :: psi_values(n_psi)  !< Poloidal flux at the surfaces
    integer,                  intent(in)        :: n_pts              !< Number of points per flux surface segment
    integer,allocatable,      intent(out)       :: n_seg(:)           !< Number of segments of each flux surface
    real*8, allocatable,      intent(out)       :: R(:,:,:)           !< R coordinates of points
    real*8, allocatable,      intent(out)       :: Z(:,:,:)           !< Z coordinates of points
    integer,                  intent(out)       :: error              !< Error flag
    
    ! --- Local variables
    integer                  :: i_elm
    integer                  :: i, j, l, n_seg_max, ip
    real*8                   :: ss1, dss1, ss2, dss2, tt1, dtt1, tt2, dtt2
    real*8                   :: u, si, dsi, ti, dti
    real*8                   :: R_s, R_t, R_st, R_ss, R_tt, Z_s, Z_t, Z_st, Z_ss, Z_tt
    type (type_surface_list) :: surface_list
    
    error = 0
    
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ! --- Find flux surfaces
    surface_list%n_psi = n_psi
    allocate (surface_list%psi_values(n_psi))
    surface_list%psi_values = psi_values
    call find_flux_surfaces(xpoint,xcase,node_list,element_list,surface_list)
    
    ! --- Determine maximum number of flux surface segments
    n_seg_max = 0
    do l = 1, surface_list%n_psi
      n_seg_max = max( n_seg_max, surface_list%flux_surfaces(l)%n_pieces )
    end do
    
    if ( allocated(R) )     deallocate(R)
    if ( allocated(Z) )     deallocate(Z)
    if ( allocated(n_seg) ) deallocate(n_seg)
    allocate (R(n_psi, n_seg_max, n_pts))
    allocate (Z(n_psi, n_seg_max, n_pts))
    !allocate (n_seg(n_psi_max)
    allocate (n_seg(n_psi))
    !allocate (n_seg(surface_list%flux_surfaces(1)%n_pieces))
    ! --- Loop over all flux surfaces
    do i=1, n_psi
      n_seg(i) = surface_list%flux_surfaces(i)%n_pieces
      
      ! --- Loop over all segments of this flux surface
      do j=1,surface_list%flux_surfaces(i)%n_pieces
        
        ! --- Bezier element, in which the current flux surface segment is located
        i_elm = surface_list%flux_surfaces(i)%elm(j)
        
        ss1  = surface_list%flux_surfaces(i)%s(1,j)
        dss1 = surface_list%flux_surfaces(i)%s(2,j)
        ss2  = surface_list%flux_surfaces(i)%s(3,j)
        dss2 = surface_list%flux_surfaces(i)%s(4,j)
        
        tt1  = surface_list%flux_surfaces(i)%t(1,j)
        dtt1 = surface_list%flux_surfaces(i)%t(2,j)
        tt2  = surface_list%flux_surfaces(i)%t(3,j)
        dtt2 = surface_list%flux_surfaces(i)%t(4,j)
        
        ! --- Loop over points in a flux surface segment
        do ip = 1, n_pts
          u = -1. + 2.*float(ip-1)/float(n_pts-1)
          
          ! --- Determine s and t values of the current point inside element i_elm
          call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
          call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
          
          ! --- Determine (R,Z)-coordinates of the current point on the current flux surface
          call interp_RZ(node_list, element_list, i_elm, si, ti, R(i, j, ip), R_s, R_t, R_st, R_ss,&
            R_tt, Z(i, j, ip), Z_s, Z_t, Z_st, Z_ss, Z_tt)
          
        end do
        
      end do
      
    end do
  
  end subroutine find_flux_surfaces2
  
  
  
  !> Implements the 'heatfluxpattern' command: Determines the heat flux pattern through an
  !! axisymmetric target plate parametrized as a set of line segments.
  recursive subroutine heatfluxpattern(command, file_handle, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    integer, parameter :: MAX_TARGET_PLATE_POINTS = 10001
    type type_target_plate
      integer :: n_points
      real*8  :: R(MAX_TARGET_PLATE_POINTS)
      real*8  :: Z(MAX_TARGET_PLATE_POINTS)
      real*8  :: length(MAX_TARGET_PLATE_POINTS)
    end type type_target_plate
    
    ! --- Local variables
    real*8, parameter   :: pi = 3.1415926535897932384626433832795029
    real*8, parameter   :: kB  = 1.3804688d-23    ! Boltzmann constant [ J / K ]
    real*8, parameter   :: mu0 = 4.d-7*pi         ! Magnetic constant [ N / A^2 ]
    character(len=1024) :: target_filename        ! File containing the target plate data
    character(len=1024) :: filename               ! Data output filename
    real*8              :: n0, rho0, m_ion        ! Core particle/mass density, ion mass
    type(type_target_plate), allocatable :: target_plate
    integer             :: ioerror
    integer             :: i, j, k
    integer             :: linepoints, tor_points, ielm, iharm
    real*8              :: phi, len, frac, xjac, four_mode
    real*8              :: R_out, Z_out, s, t
    real*8              :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    real*8              :: norm_R, norm_Z, deltaR, deltaZ
    real*8              :: dpsi_ds, dpsi_dt, du_ds, du_dt, rho, Temp, dT_ds, dT_dt
    real*8              :: dpsi_dR, dpsi_dZ, du_dR, du_dZ, dT_dR, dT_dZ
    real*8              :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8              :: dT_dR_SI, dT_dZ_SI
    real*8              :: b_R, b_Z
    real*8              :: psi, psi_norm
    real*8              :: gradpar_T_R, gradpar_T_Z
    real*8              :: velocity_R, velocity_Z
    real*8              :: heatflux_convect_R, heatflux_convect_Z, heatflux_convect_norm
    real*8              :: heatflux_diff_R, heatflux_diff_Z, heatflux_diff_norm
    real*8              :: heatflux_R, heatflux_Z, heatflux_norm
    real*8              :: ZK_prof
    real*8              :: energy_dens
    real*8              :: vpar, vpar_R_SI, vpar_Z_SI
    real*8, allocatable :: collected_data(:,:,:)
    
    error = 0
    
    ! --- Some basic checks
    call check_param_count(command, (/3,4/), error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ! --- Extract command parameters.
    target_filename = command%option(2)
    
    n0 = to_float(command%option(3), error)
    if ( error /= 0 ) return
    
    m_ion = to_float(command%option(4), error)
    if ( error /= 0 ) return
    
    rho0 = m_ion * n0
    
    ! --- Read target plate data.
    open(unit=file_handle, file=target_filename, status='old', action='read', iostat=error)
    if ( error /= 0 ) then
      write(*,*) 'ERROR in routine heatfluxpattern opening file "', trim(target_filename), '".'
      call specific_help('heatfluxpattern')
      return
    end if
    allocate( target_plate )
    target_plate%n_points = 0
    do
      read(file_handle, *, iostat=ioerror) target_plate%R(target_plate%n_points+1), &
        target_plate%Z(target_plate%n_points+1)
      if ( ioerror /= 0 ) exit
      target_plate%n_points = target_plate%n_points+1
    end do
    close(unit=file_handle)
    if ( target_plate%n_points < 2 ) then
      write(*,*) 'ERROR in routine heatfluxpattern: Could not read line segments from file "',     &
        trim(target_filename), '".'
      error = 1
      return
    end if
    
    ! --- Determine length of target plate in the poloidal cut.
    target_plate%length(:) = 0.d0
    do i = 2, target_plate%n_points
      target_plate%length(i) = target_plate%length(i-1) + sqrt( &
        (target_plate%R(i)-target_plate%R(i-1))**2 + &
        (target_plate%Z(i)-target_plate%Z(i-1))**2 )
    end do
    
    ! --- Determine heat flux onto target plate
    linepoints = get_int_setting('linepoints', error)
    tor_points = get_int_setting('tor_points', error)
    allocate( collected_data(linepoints+1, tor_points+1, 7) )
    
    do i = 1, linepoints+1 ! Points along the target in poloidal direction
      
      write(*,'(i5," of",i5)')  i, linepoints+1
      
      ! --- Determine (R,Z)-position on target and normal vector
      len = (real(i-1)/real(linepoints)) * target_plate%length(target_plate%n_points)
      do k = 2, target_plate%n_points
        if ( len <= target_plate%length(k) ) then
          ! --- (R,Z)-position on target
          frac = (len-target_plate%length(k-1)) / &
            (target_plate%length(k)-target_plate%length(k-1))
          R    = (1.d0-frac) * target_plate%R(k-1) + frac * target_plate%R(k)
          Z    = (1.d0-frac) * target_plate%Z(k-1) + frac * target_plate%Z(k)
          
          ! --- Normal vector
          deltaR = target_plate%R(k) - target_plate%R(k-1)
          deltaZ = target_plate%Z(k) - target_plate%Z(k-1)
          norm_R = + deltaZ / sqrt( deltaR**2 + deltaZ**2 )
          norm_Z = - deltaR / sqrt( deltaR**2 + deltaZ**2 )
          exit
        end if
      end do
      
      do j = 1, tor_points+1 ! Points along the target in toroidal direction
        
        ! --- Determine phi-position on target
        phi = real(j-1)/real(tor_points)
        
        ! --- Determine dR/ds, dR/dt, ..., and the 2D-Jacobian
        call find_RZ(node_list, element_list, R, Z, R_out, Z_out, ielm, s, t, error)
        if ( error /= 0 ) then
          write(*,*) 'ERROR in routine heatfluxpattern calling find_RZ. Maybe the target plate'
          write(*,*) 'is not completely inside the JOREK computational domain.'
          return
        end if
        call interp_RZ(node_list, element_list, ielm, s, t, R_out, R_s, R_t, R_st, R_ss, R_tt,    &
          Z_out, Z_s, Z_t, Z_st, Z_ss, Z_tt)
        xjac = R_s * Z_t - R_t * Z_s
        
        ! --- Determine variable values and derivatives at position (R,Z,phi).
        psi     = 0.d0
        dpsi_ds = 0.d0
        dpsi_dt = 0.d0
        du_ds   = 0.d0
        du_dt   = 0.d0
        rho     = 0.d0
        Temp    = 0.d0
        dT_ds   = 0.d0
        dT_dt   = 0.d0
        vpar    = 0.d0
        
        do iharm = 1, n_tor
          
          if ( iharm == 1 ) then
            four_mode = 1.d0
          else if ( mod(iharm,2) == 0 ) then
            four_mode = cos( mode(iharm) * phi )
          else
            four_mode = sin( mode(iharm) * phi )
          end if
          
          call interp(node_list, element_list, ielm, 1, iharm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
          psi     = psi     + P   * four_mode
          dpsi_ds = dpsi_ds + P_s * four_mode
          dpsi_dt = dpsi_dt + P_t * four_mode
          
          call interp(node_list, element_list, ielm, 2, iharm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
          du_ds   = du_ds   + P_s * four_mode
          du_dt   = du_dt   + P_t * four_mode
          
          call interp(node_list, element_list, ielm, 5, iharm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
          rho     = rho     + P   * four_mode
          
          call interp(node_list, element_list, ielm, 6, iharm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
          Temp    = Temp    + P   * four_mode
          dT_ds   = dT_ds   + P_s * four_mode
          dT_dt   = dT_dt   + P_t * four_mode
          
          if ( jorek_model >= 300 ) then
            call interp(node_list, element_list, ielm, 7, iharm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
            vpar = vpar + P * four_mode
          end if
          
        end do
        
        psi_norm = (psi - psi_axis) / (psi_bnd - psi_axis)
        dpsi_dR  = (dpsi_ds * Z_t - dpsi_dt * Z_s ) / xjac
        dpsi_dZ  = (dpsi_ds * R_t - dpsi_dt * R_s ) / xjac
        du_dR    = (du_ds   * Z_t - du_dt   * Z_s ) / xjac
        du_dZ    = (du_ds   * R_t - du_dt   * R_s ) / xjac
        dT_dR    = (dT_ds   * Z_t - dT_dt   * Z_s ) / xjac
        dT_dZ    = (dT_ds   * R_t - dT_dt   * R_s ) / xjac
        
        ! --- Temperature in SI units
        dT_dR_SI = dT_dR / ( kB * mu0 * n0 )
        dT_dZ_SI = dT_dZ / ( kB * mu0 * n0 )
        
        ! --- Magnetic field direction vector b = B / |B|
        b_R = + dpsi_dZ / sqrt( F0**2 + dpsi_dR**2 + dpsi_dZ**2 )
        b_Z = - dpsi_dR / sqrt( F0**2 + dpsi_dR**2 + dpsi_dZ**2 )
        
        ! --- Parallel velocity in SI units
        vpar_R_SI = vpar / sqrt( mu0 * rho0 ) * b_R
        vpar_Z_SI = vpar / sqrt( mu0 * rho0 ) * b_Z
        
        ! --- Parallel temperature gradient
        gradpar_T_R = b_R * ( b_R * dT_dR_SI + b_Z * dT_dZ_SI )
        gradpar_T_Z = b_Z * ( b_R * dT_dR_SI + b_Z * dT_dZ_SI )
        
        ! --- Poloidal velocity components
        velocity_R = - R * du_dZ / sqrt( mu0 * rho0 )
        velocity_Z = + R * du_dR / sqrt( mu0 * rho0 )
        
        ! --- Energy density
        energy_dens = 3.d0/2.d0 * rho * Temp / mu0
        
        ! --- Heat diffusion coefficients
        ZK_prof = ZK_perp(1) * ((1.d0-ZK_perp(2)) + ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psi_norm-ZK_perp(5))/ZK_perp(4))))
        
        ! --- Convective heat flux
        heatflux_convect_R = energy_dens * ( velocity_R + vpar_R_SI )
        heatflux_convect_Z = energy_dens * ( velocity_Z + vpar_Z_SI )
        
        ! --- Diffusive heat flux
        heatflux_diff_R = ZK_prof * dT_dR_SI + (ZK_par - ZK_prof) * gradpar_T_R
        heatflux_diff_Z = ZK_prof * dT_dZ_SI + (ZK_par - ZK_prof) * gradpar_T_Z
        
        ! --- Total heat flux
        heatflux_R = heatflux_convect_R + heatflux_diff_R
        heatflux_Z = heatflux_convect_Z + heatflux_diff_Z
        
        ! --- Heat flux through the target (normal to the target plate)
        heatflux_convect_norm = abs( heatflux_convect_R * norm_R + heatflux_convect_Z * norm_Z )
        heatflux_diff_norm    = abs( heatflux_diff_R    * norm_R + heatflux_diff_Z    * norm_Z )
        heatflux_norm         = abs( heatflux_R         * norm_R + heatflux_Z         * norm_Z )
        
        collected_data(i,j,:) = (/ R, Z, phi, len, heatflux_convect_norm, heatflux_diff_norm,      &
          heatflux_norm /)
        
      end do
      
    end do
    
    ! --- Open data output file.
    write(filename,'(a,i5.5,a)') 'heatfluxpattern_step', index_start, '.dat'
    open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
    if (error /= 0) then
      write(*,*) 'ERROR in routine heatfluxpattern: Creating/Opening file "', trim(filename), &
        '" failed.'
      return
    end if
    
    ! --- Write out data
    do j = 1, 1!tor_points+1 ! Points along the target in toroidal direction
      do i = 1, linepoints+1 ! Points along the target in poloidal direction
        write(file_handle,*) collected_data(i,j,4), collected_data(i,j,6)
        write(99,'(8es18.8)') t_start, maxval(collected_data(:,:,5)), maxval(collected_data(:,:,6)), maxval(collected_data(:,:,7)) !###
      end do
      write(file_handle,*)
      write(file_handle,*)
    end do
    
    ! --- Close data output file
    close(unit=file_handle)
    
    deallocate( collected_data )
    
  end subroutine heatfluxpattern
  
  
  
  !> Implements the 'average' command: Writes out the flux-surface averaged profile of a variable
  recursive subroutine average(command, first_step, file_handle, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    ! --- Local variables
    character(len=1024)      :: filename
    integer                  :: i_elm, ivar
    integer                  :: i, j, iopt, nplot, ip
    real*8                   :: ss1, dss1, ss2, dss2, tt1, dtt1, tt2, dtt2
    real*8                   :: u, si, dsi, ti, dti
    real*8, allocatable      :: psi_n(:)     ! Normalized poloidal flux at flux surfaces
    real*8, allocatable      :: P(:), R(:), Z(:) ! Averaged quantity and R, Z coords on pts of flux surface segment
    real*8                   :: P_s, P_t, P_st, P_ss, P_tt ! Derivatives of averaged quantity w/ respect to s and t
    real*8                   :: R_s, R_t, R_st, R_ss, R_tt ! Derivatives of R with respect to s and t
    real*8                   :: Z_s, Z_t, Z_st, Z_ss, Z_tt ! Derivatives of Z with respect to s and t
    real*8                   :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt ! Poloidal flux, Psi, and derivatives
    real*8                   :: xjac         ! 2D Jacobian
    type (type_surface_list) :: surface_list ! List of flux surfaces
    real*8                   :: dx           ! Distance between two points of a flux-surface
    real*8                   :: psi_R, psi_Z ! Derivatives of the poloidal flux w/ respect to R and Z
    real*8                   :: abs_grad_psi ! Absolute value of grad Psi required for flux surface averaging
    real*8                   :: nom          ! Nominator of the flux-surface average (oint dx P / |grad Psi|)
    real*8                   :: den          ! Denominator of the flux-surface average (oint dx 1 / |grad Psi|)
    
    error = 0
    
    ! --- Some basic checks
    call check_param_count(command, (/1,2,3,4,5,6,7,8/), error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ! --- Initialization
    nplot = 11 ! Number of points in each flux surface segment
    surface_list%n_psi = get_int_setting('surfaces', error) ! Number of flux surfaces
    allocate (surface_list%psi_values(surface_list%n_psi))
    allocate (psi_n(1:surface_list%n_psi), P(1:nplot), R(1:nplot), Z(1:nplot))
    
    ! --- Find flux surfaces at psi_n = 0.001 ... 0.999
    do i = 1, surface_list%n_psi
      psi_n(i) = 0.001d0 + 0.998d0 * real(i-1) / real(surface_list%n_psi-1)
      surface_list%psi_values(i) = psi_axis + ( psi_bnd - psi_axis ) * psi_n(i)
    end do
    call find_flux_surfaces(xpoint,xcase,node_list,element_list,surface_list)
      
    ! --- Loop over variables to average
    do iopt = 2, command%n_options
    
      ! Quantity to be averaged
      ivar = get_variable_number(command%option(iopt), error)    
      if ( error /= 0 ) return
      
      
      ! --- Open file
      write(filename,'(a,a,i5.5,a,i5.5,a)') lower_case(trim(variable_names(ivar))),                &
        '_flux-surface-averaged_steps', loop_min_step, '-', loop_max_step, '.dat'
      if ( first_step == .true. ) then
        open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
      else
        open(unit=file_handle, file=filename, status='old',     action='write', access='append',   &
          iostat=error)
      end if
      if (error /= 0) then
        write(*,*) 'ERROR in routine average: Creating/Opening file "', trim(filename), '" failed.'
        return
      end if
      write(file_handle,'(a,i5.5,a)') '# step ', index_start, ':'
      
      ! --- Loop over all flux surfaces
      do i = 1, surface_list%n_psi
        nom = 0.d0
        den = 0.d0
        
        ! --- Loop over all segments of this flux surface
        do j = 1,surface_list%flux_surfaces(i)%n_pieces
          
          ! --- Bezier element, in which the current flux surface segment is located
          i_elm = surface_list%flux_surfaces(i)%elm(j)
          
          ss1  = surface_list%flux_surfaces(i)%s(1,j)
          dss1 = surface_list%flux_surfaces(i)%s(2,j)
          ss2  = surface_list%flux_surfaces(i)%s(3,j)
          dss2 = surface_list%flux_surfaces(i)%s(4,j)
          
          tt1  = surface_list%flux_surfaces(i)%t(1,j)
          dtt1 = surface_list%flux_surfaces(i)%t(2,j)
          tt2  = surface_list%flux_surfaces(i)%t(3,j)
          dtt2 = surface_list%flux_surfaces(i)%t(4,j)
          
          ! --- Loop over nplot points within a flux surface segment
          do ip = 1, nplot
            u = -1. + 2.*float(ip-1)/float(nplot-1)
            
            ! --- Determine s and t values of the current point inside element i_elm
            call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
            call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
            
            ! --- Determine (R,Z)-coordinates of current position
            call interp_RZ(node_list, element_list, i_elm, si, ti, R(ip), R_s, R_t, R_st, R_ss, R_tt,&
              Z(ip), Z_s, Z_t, Z_st, Z_ss, Z_tt)
            
            ! --- Determine psi values and derivatives at current position
            call interp(node_list, element_list, i_elm, ivar, 1, si, ti, psi, psi_s, psi_t, psi_st,  &
              psi_ss, psi_tt)
            
            ! --- Determine averaged quantity and derivatives at current position
            call interp(node_list, element_list, i_elm, ivar, 1, si, ti, P(ip), P_s, P_t, P_st, P_ss,&
               P_tt)
             
            ! --- Determine the absolute value of grad Psi
            xjac         = R_s * Z_t - R_t * Z_s ! 2D Jacobian
            psi_R        = ( psi_s * Z_t - psi_t * Z_s ) / xjac
            psi_Z        = ( psi_s * R_t - psi_t * R_s ) / xjac
            abs_grad_psi = sqrt( psi_R**2 + psi_Z**2 )
            
            ! --- Contributions to nominator/denominator of flux surface average \<P\> defined as
            !
            !     \f$ <P> = ( oint dx P / |grad Psi| ) / ( oint dx 1 / |grad Psi| ) \f$
            !
            if ( ( ip > 1 ) .and. ( .not. in_private(R(ip),Z(ip),psi,xpoint,xcase,R_xpoint,        &
              Z_xpoint,psi_xpoint,99.d0) ) ) then
              dx  = sqrt( (R(ip)-R(ip-1))**2 + (Z(ip)-Z(ip-1))**2 )
              nom = nom + dx * ( P(ip) + P(ip-1) ) / 2.d0 / abs_grad_psi
              den = den + dx / abs_grad_psi
            end if
            
          end do ! ( loop over points within a flux surface segment)
          
        end do ! (loop over flux surface segments)
       
        ! --- Write out flux surface average at current surface to the file
        write(file_handle,*) psi_n(i), nom / den
        
      end do ! (loop over flux surfaces)
      
      ! --- Close file
      write(file_handle,*)
      write(file_handle,*)
      close (unit=file_handle)
      
    end do ! (loop over variables to average)
    
    ! --- Clean up
    deallocate (psi_n, surface_list%psi_values, P, R, Z)
    
  end subroutine average
  
  
  
  !> Implements the 'set variable value' command
  subroutine set(command, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: error   !< Error flag
      
    ! --- Local variables
    character(len=128)   :: name
    character(len=1024)  :: value
    
    error = 0
    
    ! --- Transformation of input data
    call check_param_count(command, (/0, 2/), error)
    if ( error /= 0 ) return
    
    if ( command%n_options == 1 ) then
      call print_settings()
    else
      name  = command%option(2)
      value = command%option(3)
      call set_setting(name, value, error)
    end if
    
  end subroutine set
  
  
  
  !> List all existing restart files
  subroutine timesteps()
    
    character(len=256)  :: filename
    logical             :: file_exists
    integer             :: i
    
    write(*,'(a)') 'Available restart files:'
    do i = 0, 99999
      write (filename,'(a, i5.5, a)'), 'jorek', i, '.rst'
      inquire (file=filename, exist=file_exists)
      if (file_exists) write(*,'(i6)',advance='no') i
    end do
    write(*,*)
    
  end subroutine timesteps
  
  
  
  !> Load a specific namelist input file
  subroutine load_namelist(command, file_handle, error)
    
    use phys_module     
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    ! --- Local variables
    character(len=1024) ::  filename
    logical             ::  file_exists
 
    error = 0
    
    call check_param_count(command, 1, error)
    if ( error /= 0 ) return
    
    filename = trim(command%option(2))
    inquire (file=filename, exist=file_exists)
      
    ! --- Read the input namelist file
    if (file_exists) then
      call initialise_parameters(0, filename)
      input_loaded = .true.
    else
      error = 1
      write(*,*) 'ERROR: input file "', trim(filename), '" does not exist.'
      call specific_help('namelist')
    end if

  end subroutine load_namelist
  
  
  
  !> Write out the total plasma current
  recursive subroutine global_parameters(command, first_step, file_handle, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    ! --- Local variables
    character(len=1024) :: filename
    real*8     :: Bgeo, aminor, psi_limit, current, beta_p, beta_t, beta_n
    real*8     :: density, density_in, density_out, pressure, pressure_in, pressure_out
    
    call check_param_count(command, 0, error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ! --- Open file
    write(filename,'(a,i5.5,a, i5.5, a)') 'plasma-current_steps', loop_min_step, '-',              &
      loop_max_step, '.dat'
    if ( first_step == .true. ) then
      open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
      write(file_handle,*) '# time            ###'
    else
      open(unit=file_handle, file=filename, status='old', action='write', access='append',         &
        iostat=error)
    end if
    if (error /= 0) then
      write(*,*) 'ERROR in routine global_parameters: Creating/Opening file "', trim(filename),    &
        '" failed.'
      return
    end if
    
    ! ### add area and volume as parameters to integrals and output them
    
    call integrals(node_list,element_list,R_xpoint,Z_xpoint,psi_xpoint,psi_limit,aminor,Bgeo,      &
      current,beta_p,beta_t,beta_n,density,density_in,density_out,pressure,pressure_in,pressure_out)
    
    write(file_handle,'(6ES16.7,3x,a,i5.5)') t_start, psi_limit, current/1.e6, beta_p, beta_t,     &
      beta_n, '#step', index_start
    
    close(unit=file_handle)
    
  end subroutine global_parameters
  
  
  
  !> Write out the position of the magnetic axis
  recursive subroutine magnetic_axis(command, first_step, file_handle, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag

    ! --- Local variables
    character(len=1024) :: filename
    
    error = 0
    
    call check_param_count(command, 0, error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
        
    ! --- Open file
    write(filename,'(a,i5.5,a, i5.5, a)') 'magnetic-axis_steps', loop_min_step, '-', loop_max_step,&
      '.dat'
    if ( first_step == .true. ) then
      open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
      write(file_handle,*) '# time            R_axis          Z_axis         psi_axis'
    else
      open(unit=file_handle, file=filename, status='old', action='write', access='append',         &
        iostat=error)
    end if
    if (error /= 0) then
      write(*,*) 'ERROR in routine magnetic_axis: Creating/Opening file "', trim(filename),        &
        '" failed.'
      return
    end if
    
    write(file_handle,'(4ES16.7,3x,a,i5.5)') t_start, R_axis, Z_axis, psi_axis, '#step', index_start
    
    close(unit=file_handle)
    
  end subroutine magnetic_axis
  
  
  
  !> Calculate the plasma volume (region enclosed by the plasma separatrix)
  !!
  !! The plasma cross-section is split into n_bars bars. The total plasma volume then is the sum of
  !! the bar surfaces times their respective toroidal extent: V = sum_i A_i * 2*pi*R_i
  !!
  recursive subroutine plasma_volume(command, first_step, file_handle, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    ! --- Local variables
    integer                  :: i, j, k, l, m, n
    integer                  :: n_bars
    integer,allocatable      :: n_seg(:)           ! Number of segments in each fluxsurface
    integer,allocatable      :: n_pts(:)           ! Intersections of each bar with separatrix
    real*8                   :: R_min, R_max
    real*8                   :: d_R                ! delta R between two bars
    real*8                   :: dR, dR_2, dZ, dZ_2 
    real*8                   :: Z_tmp
    real*8                   :: psi_separatrix
    real*8                   :: dV, V
    real*8                   :: dZleft, dZright, Rleft, Rright ! bars
    real*8, allocatable      :: R_b(:)     ! R over all bars
    real*8, allocatable      :: Z_b(:,:)   ! Z over all bars
    real*8, allocatable      :: R(:,:,:), Z(:,:,:)
    real*8, parameter        :: pi = 3.1415926535897932384626433832795029
    character(len=1024)      :: filename
    
    error = 0
    
    n_bars = 500 ! Number of bars to split the plasma volume into
    
    ! --- Check input data
    call check_param_count(command, 0, error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ! --- Locate separatrix
    psi_separatrix = psi_axis + (psi_bnd - psi_axis) * 0.999d0
    call find_flux_surfaces2(1, (/psi_separatrix/), 20, n_seg, R, Z, error)
    
    ! --- Locate min/max value of R and the corresponding Z value
    allocate (R_b(n_bars+1))
    allocate (n_pts(n_bars+1))
    allocate (Z_b(n_bars+1,10))
    
    R_min = 1.d99
    R_max = -1.d99
    do i = 1, n_seg(1)
      do j = 1, 20
        if ( .not. in_private(R(1,i,j),Z(1,i,j),psi_separatrix,xpoint,xcase,R_xpoint,Z_xpoint,     &
          psi_xpoint,99.d0) ) then
          if (R(1, i, j) < R_min ) then
	    R_min      = R(1, i, j)
	    R_b(1)     = R(1, i, j)
	    Z_b(1,1:2) = Z(1, i, j)
	  end if
          if (R(1, i, j) > R_max ) then
	    R_max             = R(1, i, j)
	    R_b(n_bars+1)     = R(1, i, j)
	    Z_b(n_bars+1,1:2) = Z(1, i, j)
	  end if
	end if
      end do
    end do
    n_pts(:)        = 0
    n_pts(1)        = 2
    n_pts(n_bars+1) = 2
    
    d_R = (R_max - R_min) / n_bars
    
    !--- Determine intersections of the bars with the separatrix
    do k = 2, n_bars
      R_b(k) = R_min + d_R * (k - 1)
      do l = 1, n_seg(1)
        do m = 1, 19
          if ( ( ( R(1, l, m) <= R_b(k) ) .and. ( R_b(k) < R(1, l, m+1) ) ) .or.  &
	       ( ( R(1, l, m) >= R_b(k) ) .and. ( R_b(k) > R(1, l, m+1) ) ) ) then 
	    dR    = R(1, l, m+1) - R(1, l, m) ! R-distance between to separatrix points
	    dR_2  = R(1, l, m+1) - R_b(k)     ! R-distance between separatrix point and bar
	    dZ    = Z(1, l, m+1) - Z(1, l, m) ! Z-distance between separatrix points
	    dZ_2  = (dR_2 / dR) * dZ
            Z_tmp = Z(1, l, m+1) - dZ_2       ! Z-position of separatrix-bar intersection
            
	    if ( in_private(R(1,l,m),Z_tmp,psi_separatrix,xpoint,xcase,R_xpoint,Z_xpoint,          &
              psi_xpoint,99.d0) ) cycle  ! discard points on divertor legs
	    
	    n_pts(k) = n_pts(k) + 1
	    Z_b(k, n_pts(k)) = Z_tmp
            
            write(38,'(2ES16.7, a)') R_b(k), Z_tmp
          endif
        enddo
      enddo
    enddo
    
    ! --- Try to reduce number of cut-points of each bar with the separatrix to two if necessary
    !     (more than two intersections can sometimes occur if the separatrix is slightly deformed)
    do k = 1, n_bars+1
      if (n_pts(k) > 2) then
        do i = 1, 2
	  j = i + 1
          do
	    if ( j > n_pts(k) ) exit
            if ( abs( Z_b(k, j) - Z_b(k, i) ) < 0.005 ) then
	      Z_b(k,j:n_pts(k)-1) = Z_b(k,j+1:n_pts(k)) 
	      n_pts(k) = n_pts(k) - 1
	      write(*,'(a,i5)') 'WARNING: Removing redundant intersection point in bar', k
	    else
	      j = j + 1
	    end if
          end do
        end do
        if ( n_pts(k) > 2 ) then
          write(*,*) 'ERROR in routine plasma_volume: Intersection points /= 2 for bar ',k
	  error = 1
	  return
	end if
      end if
    end do
    
    ! --- Calculate the plasma volume enclosed by the separatrix
    V  = 0
    do n = 1, n_bars
      dZleft  = abs( Z_b(n, 1) - Z_b(n, 2) )
      dZright = abs( Z_b(n+1, 1) - Z_b(n+1, 2) )
      Rleft   = R_b(n)
      Rright  = R_b(n+1)
      
      !dV= average-bar-height    * bar-width      * average-toroidal-length
      dV = (dZright+dZleft)/2.d0 * (Rright-Rleft) * 2*pi*(Rright+Rleft)/2.d0
      V  = V + dV
    end do
    write(*,'(a,f8.3,a)') 'Plasma volume: ', V, ' m^3'
    
    ! --- Create plasma-volume data-file
    write(filename,'(a,i5.5,a, i5.5, a)') 'plasma-volume_steps', loop_min_step, '-', loop_max_step,&
       '.dat'
    if ( first_step == .true. ) then
      open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
    else
      open(unit=file_handle, file=filename, status='old', action='write', access='append',         &
        iostat=error)
    end if
    if (error /= 0) then
      write(*,*) 'ERROR in routine plasma_volume creating/opening file "', trim(filename), '".'
      return
    end if
    write(file_handle,'(2ES16.7 )') t_start, V
    close(unit=file_handle)
    
    deallocate (R_b, Z_b, n_pts, n_seg, R, Z)
    
  end subroutine plasma_volume
  
  
  
  !> Variable-value at a certain position
  recursive subroutine point(command, first_step, file_handle, error)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag
    
    ! --- Local variables
    integer :: ivar
    real*8  :: R
    real*8  :: Z
    real*8  :: phi
    character(len=1024) :: filename
    
    error = 0
    
    ! --- Transformation of input data
    call check_param_count(command, 4, error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    ivar = get_variable_number(command%option(2), error)
    if ( error /= 0 ) return
    
    R = to_float(command%option(3), error)
    if ( error /= 0 ) return
    
    Z = to_float(command%option(4), error)
    if ( error /= 0 ) return
    
    phi = to_float(command%option(5), error)
    if ( error /= 0 ) return
    
    if ( (R<0.d0) .or. (phi<0.d0) ) then
      write(*,*) 'ERROR in routine point: R and phi must be positive.'
      error = 1
      return
    end if
    
    100 format(a,a,sp,es10.3,a,es10.3,a,es10.3,a,ss,i5.5,a,i5.5,a)
    write(filename,100) lower_case(trim(variable_names(ivar))), '_at_', R, ',', Z, ',', phi,        &
      '_steps', loop_min_step, '-', loop_max_step, '.dat'
    
    if ( first_step == .true. ) then
      open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
    else
      open(unit=file_handle, file=filename, status='old', action='write', access='append',         &
        iostat=error)
    end if
    
    if (error /= 0) then
      write(*,*) 'ERROR in routine point: Creating/Opening file "', trim(filename), '" failed.'
      return
    end if
    write(file_handle,'(2ES16.7)') t_start, variable_value(ivar, R, Z, phi, error)
    if (error /= 0) write(*,*) 'ERROR executing variable_value.'
    
    close(unit=file_handle)
    
  end subroutine point
  
  
  
  !> Variable value along a straight line in (R,Z,phi)
  recursive subroutine line(command, first_step, file_handle, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: error       !< Error flag

    
    ! --- Local variables
    integer :: ivar, i
    character(len=128)  :: var(2)
    character(len=1024) :: filename
    real*8  :: R0, Z0, phi0    ! Start point of line
    real*8  :: R1, Z1, phi1    ! End point of line
    real*8  :: R_i, Z_i, phi_i ! Current point
    integer :: n_points  ! Number of points
    integer :: i_point   ! Current point
    real*8  :: s_point   ! Normalized length along the line (0...1)
    real*8  :: val(2)
    
    error = 0
    
    ! --- Transformation of input data
    call check_param_count(command, 8, error)
    if ( error /= 0 ) return
    call check_step_imported(error)
    if ( error /= 0 ) return
    
    do i = 1, 2
      var(i) = lower_case(command%option(i+1))
      
      if ( (var(i) /= 'r') .and. (var(i) /= 'z') .and. (var(i) /= 'phi') .and.                     &
        (var(i) /= 'normalized_length') .and. (var(i) /= 'psi_n') ) then
        
        ivar = get_variable_number(command%option(i+1), error)
        if ( error /= 0 ) return
        var(i) = variable_names(ivar)
        
      end if
      
    end do
    
    R0 = to_float(command%option(4), error)
    if ( error /= 0 ) return
    Z0 = to_float(command%option(5), error)
    if ( error /= 0 ) return
    phi0 = to_float(command%option(6), error)
    if ( error /= 0 ) return
    
    R1 = to_float(command%option(7), error)
    if ( error /= 0 ) return
    Z1 = to_float(command%option(8), error)
    if ( error /= 0 ) return
    phi1 = to_float(command%option(9), error)
    if ( error /= 0 ) return
    
    if ( (R0<0.d0) .or. (R1<0.d0) .or. (phi0<0.d0) .or. (phi1<0.d0) ) then
      write(*,*) 'ERROR in routine line: R and phi must be positive.'
      error = 1
      return
    end if
    
    ! ---- File creation
    110 format(4a,sp,6(es10.3,a),ss,2(i5.5,a))
    write(filename,110) trim(lower_case(var(1))), '_versus_', trim(var(2)), '_from_', R0, ',', Z0, &
      ',', phi0, '_to_', R1, ',', Z1, ',', phi1, '_steps', loop_min_step, '-', loop_max_step, '.dat'
    
    if ( first_step == .true. ) then !> default: true
      open(unit=file_handle, file=filename, status='replace', action='write', iostat=error)
    else
      open(unit=file_handle, file=filename, status='old',     action='write', access='append',     &
        iostat=error)
    end if
    write(file_handle,'(a,i5.5,a)') '# step ', index_start, ':'
    
    if (error /= 0) then
      write(*,*) 'ERROR in routine line: Creating/Opening file', trim(filename), '" failed.'
      return
    end if
    
    ! --- Write data to file
    n_points = get_int_setting('linepoints', error) ! get from settings: linepoints
    
    do i_point = 0, n_points
      s_point = real(i_point) / real(n_points)
      R_i = R0 + (R1-R0)*s_point
      Z_i = Z0 + (Z1-Z0)*s_point
      phi_i = phi0 + (phi1-phi0)*s_point
      do i = 1, 2
        select case ( var(i) )
          case( 'r' )
            val(i) = R_i
          case( 'z' )
            val(i) = Z_i
          case( 'phi' )
            val(i) = phi_i
          case( 'normalized_length' )
            val(i) = s_point
          case( 'psi_n' )
            val(i) = (variable_value(1, R_i, Z_i, phi_i, error) - psi_axis) / (psi_bnd - psi_axis)
          case default
            val(i) = variable_value(get_variable_number(var(i), error), R_i, Z_i, phi_i, error)
        end select
      end do
      write(file_handle,'(2ES15.7)') val(2), val(1)
      if (error /= 0) then
        write(*,*) 'ERROR executing variable_value.'
        exit
      end if
    end do
    
    write(file_handle,*)
    write(file_handle,*)
    
    close (unit=file_handle)
    
  end subroutine line
  
  
  
  !> Check if a restart file has already been imported
  recursive subroutine check_step_imported(error)
    integer, intent(out) :: error !< Error flag
    
    error = 0
  
    if ( .not. step_imported ) then
      error = 1
      write(*,*) 'ERROR: No restart file has been imported yet. Use the "for" loop:'
      call specific_help('for')
    end if
    
  end subroutine check_step_imported
  
  
  
  !> Print usage information
  subroutine help(command, error)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command !< Command to be executed
    integer,            intent(out) :: error   !< Error flag
    
    ! --- Local variables
    integer :: i
    
    error = 0
    
    if ( command%n_options == 1 ) then
      call general_help()
    else
      do i = 2, command%n_options
        call specific_help(command%option(i))
      end do
    end if
    
  end subroutine help
  
  
  
end module exec_commands

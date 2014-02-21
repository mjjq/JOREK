!> Module for execution of user commands (part of jorek2_postproc)
module exec_commands
  
  use constants
  use parameters
  use phys_module
  use data_structure
  use equil_info
  use nodes_elements
  use boundary
  use basis_at_gaussian
  use mod_new_diag
  use parse_commands
  use settings
  use convert_character
  use postproc_help
  use domains
  
  
  
  implicit none
  
  
  
  integer, parameter :: NORMAL_MODE = 1 !< Normal mode
  integer, parameter :: LOOP_MODE   = 2 !< Mode started by 'for' and ended by 'done' commands
  integer :: exec_mode = NORMAL_MODE    !< Current operation mode (NORMAL_MODE or LOOP_MODE)
  integer :: loop_min_step              !< Smallest timestep index for current loop
  integer :: loop_max_step              !< Largest timestep index for current loop
  
  integer, parameter :: MAX_QUEUE_LENGTH  = 9999         !< Maximum length of command queue
  integer            :: n_queued_commands = 0            !< Number of commands in the queue
  type(type_command) :: command_queue(MAX_QUEUE_LENGTH)  !< Queued commands
  
  logical,             private, save :: input_loaded  = .false. !< Has an input file been loaded?
  logical,             private, save :: step_imported = .false. !< Has a restart file been imported?
  type(t_equil_state), private, save :: eq !< Equilibrium state; updated when time step is loaded
  type(t_expr_list),   private, save :: expr_list
  real*8, allocatable, private, save :: result(:,:,:,:), res2d(:,:,:), res1d(:,:), res0d(:)
  
  
  
  private
  public exec_command, general_help, specific_help
  
  
  
  save
  
  
  
  contains
  
  
  
  !> Execute a command
  subroutine exec_command(command, first_step, file_handle, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    if ( get_setting('debug',ierr) == 'true' ) then
      write(*,'(a)') 'Exec_command was called with:'
      call print_command(command)
      write(*,*)
    end if
    
    ! --- In normal mode, some commands are directly executed
    if ( exec_mode == NORMAL_MODE ) then
      
      if ( ( .not. input_loaded ) .and. ( trim(command%args(0)) == 'for' ) ) then
        write(*,*) 'ERROR: No namelist input file loaded.'
        call specific_help('namelist')
        ierr = 1
        return
      end if
      
      write(*,*)
      write(*,*) '- Executing "'//trim(command%args(0))//'"'
      
      select case ( trim(command%args(0)) )
        case ( 'average' )
          call average(command, first_step, ierr)
        case ( 'global_parameters' )
!          call global_parameters(command, first_step, file_handle, ierr)
        case ( 'expressions' )
          call expressions(command, ierr)
        case ( 'fluxsurfaces' )
!          call fluxsurfaces(command, file_handle, ierr)
        case ( 'for' )
          call loop_start(command, ierr)
        case ( 'gourdon' )
!          call gourdon(command, file_handle, ierr)
        case ( 'help' )
          call help(command, ierr)
        case ( 'line' )
!          call line(command, first_step, file_handle, ierr)
        case ( 'mark_coords' )
          call mark_coords(command, ierr)
        case ( 'midplane' )
          call midplane(command, first_step, ierr)
        case ( 'namelist' )
          call load_namelist(command, file_handle, ierr)
        case ( 'params' )
          call log_parameters(0)
        case ( 'point' )
          call point(command, first_step, ierr)
        case ( 'qprofile' )
!          call qprofile(command, first_step, file_handle, ierr)
        case ( 'set' )
          call set(command, ierr)
        case ( 'timesteps' )
          call timesteps 
        case default
          write(*,*) 'Command "', trim(command%args(0)), '" does not exist'
          call general_help() 
      end select
      
    ! --- In loop mode, commands are queued and afterwards executed for each timestep separately
    else if ( exec_mode == LOOP_MODE ) then
      
      select case ( trim(command%args(0)) )
        case ( 'expressions', 'mark_coords', 'global_parameters', 'midplane', 'average', 'point' )
          call add_to_command_queue(command, ierr)
        case ( 'help' )
          call help(command, ierr)
        case ( 'done' )
          call loop_end(command, ierr)
        case default
          write(*,*) 'Command "', trim(command%args(0)), '" unknown or invalid inside a loop'
          call general_help() 
      end select
        
    end if
    
  end subroutine exec_command
  
  
  
  !> Loads a time step from a restart file if the restart file exists
  subroutine load_step(istep, ierr)
    
    ! --- Routine parameters
    integer,            intent(in)  :: istep !< Load this time step
    integer,            intent(out) :: ierr !< Error flag
    
    character(len=64) :: file_name
    logical           :: file_exists
    integer           :: ifail
    
    ierr = 0
    
    write(file_name,'(a,i5.5,a)') 'jorek', istep, '.rst'
    
    inquire(file=file_name, exist=file_exists)
    if ( .not. file_exists ) then
      ierr = 1
      return
    end if
    
    write(*,*)
    write(*,'(a,i5.5,a)') '#################### TIME STEP ', istep, ' ####################'
    write(*,*)
    
    ! --- Load the restart file
    call import_restart(node_list, element_list, trim(file_name), rst_format, ierr)
    if ( ierr /= 0 ) return
    call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
    
    ! --- Locate magnetic axis and X-point.
    call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, eq)
    call print_equil_state(eq, .false.)
    
    step_imported = .true.
    
  end subroutine load_step
  
  
  
  !> Is called at the start of a for loop.
  !!
  !! It switches the module behaviour (exec_mode) to LOOP_MODE, i.e., all commands inside the for
  !! loop are not executed immediately but collected in a command queue. All commands of the
  !! command queue are executed for each time step after the 'done' command of the for loop
  !! has been entered.
  subroutine loop_start(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ierr = 0
    
    !### check param count
    
    if ( trim(command%args(1)) /= 'step' ) then
      write(*,*) 'ERROR: Wrong syntax for "for" statement.'
      call specific_help('for')
      return
    end if
    
    loop_min_step = to_int(command%args(2),ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR: Wrong syntax for "for" statement.'
      call specific_help('for')
      return
    end if
    
    if ( command%n_args == 3 ) then ! for step <XXX> do
      if ( trim(command%args(3)) /= 'do' ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
      loop_max_step = loop_min_step
    else ! for step <XXX> to <YYY> do
      if ( trim(command%args(3)) /= 'to' ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
      loop_max_step = to_int(command%args(4),ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
      if ( trim(command%args(5)) /= 'do' ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if
    end if
    
    exec_mode = LOOP_MODE
    
  end subroutine loop_start
  
  
  
  !> Is called upon the 'done' command of a for loop and executes all commands enclosed in the loop.
  subroutine loop_end(command, ierr)
    
    include 'omp_lib.h'
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ! --- Local variables
    integer :: jcmd, istep, load_error
    integer :: thread_id, file_handle
    logical :: first_step ! Is true for the first timestep loaded in the for-loop
    
    ierr = 0
    exec_mode = NORMAL_MODE    
    
    if ( n_queued_commands == 0 ) then
      write(*,*) 'WARNING: No commands in for-loop.'
      return
    end if
    
    first_step = .true.
    do istep = loop_min_step, loop_max_step
      call load_step(istep, load_error)
      if ( load_error /= 0 ) cycle
      
      do jcmd = 1, n_queued_commands
#ifdef _OPENMP        
        thread_id   = omp_get_thread_num()
#else
        thread_id   = 0
#endif
        file_handle = 17 + thread_id
        
        call exec_command(command_queue(jcmd), first_step, file_handle, ierr)  
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR executing the following command (ignoring it):'
          call print_command(command_queue(jcmd))
          ierr = 0
        end if
        
      end do
      
      first_step = .false.
    end do
    
    if ( first_step ) then
      write(*,'(a,i5.5,a,i5.5,a)') 'WARNING: There were no restart files for steps ',              &
        loop_min_step, ' to ', loop_max_step, '.'
    end if
    
    n_queued_commands = 0
    
  end subroutine loop_end
  
  
  
  !> Add a command to the command queue (for loop)
  subroutine add_to_command_queue(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ierr = 0
    
    if ( n_queued_commands + 1 > MAX_QUEUE_LENGTH) then
      write(*,*) 'ERROR: Too many commands in the command queue of the for loop.'
      ierr = 1
      return
    end if
    
    n_queued_commands = n_queued_commands + 1
    command_queue(n_queued_commands) = command
    
  end subroutine add_to_command_queue
  
  
  
  !> Retrieve a command from the command queue (currently not used!)
  subroutine get_from_command_queue(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(out)    :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ierr = 0
    
    if ( n_queued_commands < 1 ) then
      write(*,*) 'ERROR: Cannot get a command from an empty queue.'
      ierr = 1
      return
    end if
    
    n_queued_commands = n_queued_commands - 1
    command = command_queue(1)
    command_queue(1:n_queued_commands) = command_queue(2:n_queued_commands+1)
    
  end subroutine get_from_command_queue
  
  
  
  !> Implements the 'set variable value' command
  subroutine set(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
      
    ! --- Local variables
    character(len=128)   :: name
    character(len=1024)  :: value
    
    ierr = 0
    
    ! --- Transformation of input data
    !### check param count
    
    if ( command%n_args == 0 ) then
      call print_settings()
    else
      name  = command%args(1)
      value = command%args(2)
      call set_setting(name, value, ierr)
    end if
    
  end subroutine set
  
  
  
  !> List all existing restart files
  subroutine timesteps()
    
    character(len=256)  :: filename
    logical             :: file_exists
    integer             :: i
    
    write(*,'(a)') 'Available restart files:'
    do i = 0, 99999
      write (filename,'(a, i5.5, a)') 'jorek', i, '.rst'
      inquire (file=filename, exist=file_exists)
      if (file_exists) write(*,'(i6)',advance='no') i
    end do
    write(*,*)
    
  end subroutine timesteps
  
  
  
  !> Load a specific namelist input file
  subroutine load_namelist(command, file_handle, ierr)
    
    use phys_module     
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(in)  :: file_handle !< File handle
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    character(len=1024) ::  filename
    logical             ::  file_exists
 
    ierr = 0
    
    !### check param count
    
    filename = trim(command%args(1))
    inquire (file=filename, exist=file_exists)
      
    ! --- Read the input namelist file
    if (file_exists) then
      call initialise_parameters(0, filename)
      input_loaded = .true.
    else
      ierr = 1
      write(*,*) 'ERROR: input file "', trim(filename), '" does not exist.'
      call specific_help('namelist')
    end if

  end subroutine load_namelist
  
  
  
  !> Check if a restart file has already been imported
  subroutine check_step_imported(ierr)
    integer, intent(out) :: ierr !< Error flag
    
    ierr = 0
  
    if ( .not. step_imported ) then
      ierr = 1
      write(*,*) 'ERROR: No restart file has been imported yet. Use the "for" loop:'
      call specific_help('for')
    end if
    
  end subroutine check_step_imported
  
  
  
  !> Print usage information
  subroutine help(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command !< Command to be executed
    integer,            intent(out) :: ierr    !< Error flag
    
    ! --- Local variables
    integer :: i
    
    ierr = 0
    
    if ( command%n_args == 0 ) then
      call general_help()
    else
      do i = 1, command%n_args
        call specific_help(command%args(i))
      end do
    end if
    
  end subroutine help
  
  
  
  !> Report an error message.
  subroutine report_error(routine, message, command)
    
    ! --- Routine parameters.
    character(len=*),             intent(in) :: routine
    character(len=*),             intent(in) :: message
    type(type_command), optional, intent(in) :: command
    
    write(*,*)
    write(*,*) 'ERROR in '//trim(routine)//': '//trim(message)
    write(*,*)
    
    if ( present(command) ) call specific_help(trim(command%args(0)))
    
  end subroutine report_error
  
  
  
  !> List or Select Available Expressions.
  subroutine expressions(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    character(len=1024) :: filename
    
    ierr = 0
    
    if ( command%n_args == 0 ) then
      
      call print_exprs(exprs_all)
      
    else
      
      expr_list = exprs(command%args(1:command%n_args), command%n_args)
      
      call print_exprs(expr_list,.true.)
      
    end if
    
  end subroutine expressions
  
  
  
  !> Mark some expressions as coordinates.
  subroutine mark_coords(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: n_coord
    character(len=1024) :: filename
    
    ierr = 0
    
    if ( command%n_args /= 1 ) then
      call report_error('mark_coords', 'Wrong number of parameters.', command)
      ierr = 1
      return
    end if
    
    n_coord = to_int(command%args(1), ierr)
    expr_list%n_coord = n_coord
    
  end subroutine mark_coords
  
  
  
  !> Auxilliary routine for file name construction.
  character(len=64) function step_range_string(min_step, max_step)
    
    integer, intent(in) :: min_step, max_step
    
    if ( loop_min_step /= loop_max_step ) then
      write(step_range_string,'(a,i5.5,a,i5.5)') '_steps', min_step, '-', max_step
    else
      write(step_range_string,'(a,i5.5)') '_step', min_step
    end if
    
  end function step_range_string
  
  
  
  !> Evaluate expressions at a single point.
  subroutine point(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: R, Z, phi
    integer :: units
    character(len=1024) :: filename
    
    ierr = 0
    
    if ( command%n_args /= 3 ) then
      call report_error('point', 'Wrong number of parameters.', command)
      ierr = 1
      return
    end if
    
    R = to_float(command%args(1), ierr)
    if ( ierr /= 0 ) return
    
    Z = to_float(command%args(2), ierr)
    if ( ierr /= 0 ) return
    
    phi = to_float(command%args(3), ierr)
    if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    
    100 format(a,sp,es10.3,a,es10.3,a,es10.3,2a)
    write(filename,100) 'exprs_at_R', R, '_Z', Z, '_phi', phi,                                     &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    call eval_expr(eq, units, expr_list, pol_pos(node_list,element_list,eq,R=R,Z=Z),               &
      tor_pos(phi=phi), result, ierr)
    
    call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
    
    call write_ascii_0d(ierr, eq, expr_list, res0d, FORM_TABLE, header=first_step,                 &
      filename=trim(filename), append=(.not.first_step), blanks=.false.)
    
  end subroutine point
  
  
  
  !> Evaluate expressions on the midplane.
  subroutine midplane(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: R, Z, phi
    integer :: units, npts
    character(len=1024) :: filename
    
    ierr = 0
    
    if ( command%n_args /= 0 ) then
      call report_error('midplane', 'Wrong number of parameters.', command)
      ierr = 1
      return
    end if
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('linepoints', ierr)
    
    write(filename,'(3a)') 'exprs_midplane', trim(step_range_string(loop_min_step,loop_max_step)), &
      '.dat'
    
    call midplane_profile(node_list, element_list, eq, units, expr_list, res1d, BOTH_SIDES, npts,  &
      ierr, filename=trim(filename), append=(.not.first_step) )
    
  end subroutine midplane
  
  
  
  !> Toroidally and poloidally averaged expressions.
  subroutine average(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: R, Z, phi
    integer :: units, npts
    character(len=1024) :: filename
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    
    ierr = 0
    
    if ( command%n_args /= 0 ) then
      call report_error('average', 'Wrong number of parameters.', command)
      ierr = 1
      return
    end if
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename,'(3a)') 'exprs_averaged',                                                       &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    pol_pos_list = pol_pos(node_list, element_list, eq, nPsiN=npts, nTht=6*4*n_plane) !###
    tor_pos_list = tor_pos(nphi=4*n_plane) !###
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call apply_four_filter(result, simple_filter(m=0,n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                     &
      filename=trim(filename), append=(.not.first_step), blanks=.true.)
    
  end subroutine average
  
  
  
end module exec_commands

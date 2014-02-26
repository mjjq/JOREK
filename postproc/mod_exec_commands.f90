!> Module for execution of user commands (used by jorek2_postproc)
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
  use domains
  use parse_commands
  use settings
  use convert_character
  use postproc_help
  
  
  
  
  
  implicit none
  
  
  
  
  character(len=11), parameter, private :: DIR = './postproc/' !< Output goes into this directory!
  
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
  logical,             private, save :: dir_created   = .false. !< Postproc directory created?
  logical,             private, save :: verbose
  logical,             private, save :: debug
  type(t_equil_state), private, save :: eq !< Equilibrium state; updated when time step is loaded
  type(t_expr_list),   private, save :: expr_list
  real*8, allocatable, private, save :: result(:,:,:,:), res2d(:,:,:), res1d(:,:), res0d(:)
  complex*16, allocatable, private, save :: cp(:,:,:,:)
  
  
  
  
  
  private
  public exec_command, general_help, specific_help, clean_up
  
  
  
  
  
  save
  
  
  
  
  
  contains
  
  
  
  
  
  !> Execute a command
  subroutine exec_command(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    if ( .not. dir_created ) then
      call system('mkdir -p '//DIR)
      dir_created = .true.
    end if
    
    if ( get_setting('debug',ierr) == 'true' ) then
      write(*,'(a)') 'Exec_command was called with:'
      call print_command(command)
      write(*,*)
    end if
    
    ! --- Free shared module arrays
    call clean_up()
    
    ! --- In normal mode, some commands are directly executed
    if ( exec_mode == NORMAL_MODE ) then
      
      if ( ( .not. input_loaded ) .and. ( trim(command%args(0)) == 'for' ) ) then
        write(*,*) 'ERROR: No namelist input file loaded.'
        call specific_help('namelist')
        ierr = 1
        return
      end if
      
      verbose = get_log_setting('verbose', ierr)
      debug   = get_log_setting('debug', ierr)
      
      if ( verbose .or. debug ) then
        write(*,*)
        write(*,*) '- Executing "'//trim(command%args(0))//'"'
      end if
      
      select case ( trim(command%args(0)) )
        case ( 'average' )
          call average(command, first_step, ierr)
!        case ( 'global_params' )
!          call global_params(command, first_step, ierr)
        case ( 'equil_params' )
          call equil_params(command, first_step, ierr)
        case ( 'expressions' )
          call expressions(command, ierr)
        case ( 'fluxsurfaces' )
          call fluxsurfaces(command, ierr)
        case ( 'for' )
          call loop_start(command, ierr)
        case ( 'four2d' )
          call four2d(command, ierr)
        case ( 'gourdon' )
          call gourdon(command, first_step, ierr)
        case ( 'help' )
          call help(command, ierr)
        case ( 'jorek-units' )
          call select_jorek_units(command, ierr)
        case ( 'pol_line' )
          call pol_line(command, first_step, ierr)
        case ( 'tor_line' )
          call tor_line(command, first_step, ierr)
        case ( 'mark_coords' )
          call mark_coords(command, ierr)
        case ( 'midplane' )
          call midplane(command, first_step, ierr)
        case ( 'namelist' )
          call load_namelist(command, ierr)
        case ( 'params' )
          call log_parameters(0)
        case ( 'point' )
          call point(command, first_step, ierr)
        case ( 'qprofile' )
          call qprofile(command, first_step, ierr)
        case ( 'separatrix' )
          call separatrix(command, ierr)
        case ( 'set' )
          call set(command, ierr)
        case ( 'si-units' )
          call select_si_units(command, ierr)
        case ( 'timesteps' )
          call timesteps 
        case default
          write(*,*) 'Command "', trim(command%args(0)), '" does not exist'
          call general_help() 
      end select
      
    ! --- In loop mode, commands are queued and afterwards executed for each timestep separately
    else if ( exec_mode == LOOP_MODE ) then
      
      select case ( trim(command%args(0)) )
        case ( 'expressions', 'mark_coords', 'global_parameters', 'midplane', 'average', 'point',  &
          'pol_line', 'tor_line', 'equil_params', 'qprofile', 'fluxsurfaces', 'separatrix',        &
          'four2d', 'gourdon', 'jorek-units', 'si-units' )
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
    
    t_now         = t_start
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
    
    ! --- Some checks.
    call check_args(command%n_args,ierr,3,5);  if ( ierr /= 0 ) return
    
    if ( trim(command%args(1)) /= 'step' ) then
      call report_error('for', 'Wrong syntax.', command)
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
        
        call exec_command(command_queue(jcmd), first_step, ierr)  
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR executing the following command (ignoring it):'
          call print_command(command_queue(jcmd))
          call specific_help(command_queue(jcmd)%args(0))
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
  
  
  
  
  
  !> Clean up global module arrays.
  subroutine clean_up()
    if ( allocated(result) ) deallocate(result)
    if ( allocated(res2d)  ) deallocate(res2d)
    if ( allocated(res1d)  ) deallocate(res1d)
    if ( allocated(res0d)  ) deallocate(res0d)
    if ( allocated(cp)     ) deallocate(cp)
  end subroutine clean_up
  
  
  
  
  
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
  subroutine load_namelist(command, ierr)
    
    use phys_module     
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    character(len=1024) ::  filename
    logical             ::  file_exists
    
    ierr = 0
    
    ! --- Some checks.
    call check_args(command%n_args,ierr,1);  if ( ierr /= 0 ) return
    
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
  
  
  
  
  
  !> Check if expressions have been selected.
  subroutine check_exprs_selected(ierr)
    integer, intent(out) :: ierr !< Error flag
    
    ierr = 0
    if ( expr_list%n_expr <= 0 ) then
      ierr = 1
      write(*,*) 'ERROR: No physical expressions selected yet. Use command "expressions":'
      call specific_help('expressions')
    end if
    
  end subroutine check_exprs_selected
  
  
  
  
  
  !> Check number of arguments correct.
  subroutine check_args(nargs,ierr,ok1,ok2,ok3,ok4)
    integer, intent(in)              :: nargs               !< Number of arguments.
    integer, intent(inout)           :: ierr                !< Error code.
    integer, intent(in),    optional :: ok1, ok2, ok3, ok4
    
    ierr = 1
    if ( present(ok1) .and. (nargs==ok1) ) ierr = 0
    if ( present(ok2) .and. (nargs==ok2) ) ierr = 0
    if ( present(ok3) .and. (nargs==ok3) ) ierr = 0
    if ( present(ok4) .and. (nargs==ok4) ) ierr = 0
    if ( ierr /= 0 ) write(*,*) 'ERROR: Wrong number of parameters for command.'
  end subroutine check_args
  
  
  
  
  
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
    
    ierr = 0
    
    call check_args(command%n_args,ierr,1);  if ( ierr /= 0 ) return
    n_coord = to_int(command%args(1), ierr); if ( ierr /= 0 ) return
    expr_list%n_coord = n_coord
    
  end subroutine mark_coords
  
  
  
  
  
  !> Auxilliary routine for file name construction.
  character(len=64) function step_range_string(min_step, max_step)
    
    integer, intent(in) :: min_step, max_step
    
    if ( loop_min_step /= loop_max_step ) then
      write(step_range_string,'(a,i5.5,a,i5.5)') '_s', min_step, '..', max_step
    else
      write(step_range_string,'(a,i5.5)') '_s', min_step
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
    
    ! --- Some checks
    call check_args(command%n_args,ierr,3);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    R     = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Z     = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    phi   = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
    units = get_int_setting('units', ierr)
    
    write(filename,'(9a)') DIR, 'exprs_at_R', trim(real2str(R)), '_Z', trim(real2str(Z)), '_p',    &
      trim(real2str(phi)), trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
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
    integer :: units, npts
    character(len=1024) :: filename, comment
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('linepoints', ierr)
    
    write(filename,'(4a)') DIR, 'exprs_midplane',                                                  &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    write(comment,'(a,i6.6)') 'time step #', index_start
    
    call midplane_profile(node_list, element_list, eq, units, expr_list, res1d, BOTH_SIDES, npts,  &
      ierr, filename=trim(filename), append=(.not.first_step), comment=trim(comment) )
    
  end subroutine midplane
  
  
  
  
  
  !> Expressions along a line in the poloidal plane.
  subroutine pol_line(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: Rstart, Zstart, Rend, Zend, phi
    integer :: units, npts
    character(len=1024) :: filename, comment
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,5);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    Rstart = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Zstart = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    Rend   = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
    Zend   = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    phi    = to_float(command%args(5), ierr); if ( ierr /= 0 ) return
    units  = get_int_setting('units', ierr);      if ( ierr /= 0 ) return
    npts   = get_int_setting('linepoints', ierr); if ( ierr /= 0 ) return
    
    write(filename,'(15a)') DIR, 'exprs_along_line_R', trim(real2str(Rstart)), '..',               &
      trim(real2str(Rend)), '_Z', trim(real2str(Zstart)), '..', trim(real2str(Zend)), '_p',        &
      trim(real2str(phi)), trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    write(comment,'(a,i6.6)') 'time step #', index_start
    
    call pol_lineout(node_list, element_list, eq, units, expr_list, res1d, phi, Rstart, Zstart,    &
      Rend, Zend, npts, ierr, filename, append=(.not.first_step), comment=trim(comment) )
    
  end subroutine pol_line
  
  
  
  
  
  !> Expressions along a toroidal line.
  subroutine tor_line(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: R, Z, phi_start, phi_end
    integer :: units, npts
    character(len=1024) :: filename, comment
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,4);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    R         = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Z         = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    phi_start = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
    phi_end   = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    units     = get_int_setting('units', ierr);      if ( ierr /= 0 ) return
    npts      = get_int_setting('linepoints', ierr); if ( ierr /= 0 ) return
    
    write(filename,'(15a)') DIR, 'exprs_along_line_R', trim(real2str(R)), '_Z', trim(real2str(Z)), &
      '_p', trim(real2str(phi_start)), '..', trim(real2str(phi_end)),                              &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    write(comment,'(a,i6.6)') 'time step #', index_start
    
    call tor_lineout(node_list, element_list, eq, units, expr_list, res1d, phi_start, phi_end, R,  &
      Z, npts, ierr, filename, append=(.not.first_step), comment=trim(comment) )
    
  end subroutine tor_line
  
  
  
  
  
  !> Toroidally and poloidally averaged expressions.
  subroutine average(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts
    character(len=1024) :: filename, comment
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename,'(4a)') DIR, 'exprs_averaged',                                                  &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    pol_pos_list = pol_pos(node_list, element_list, eq, nPsiN=npts, nTht=6*4*n_plane) !###
    tor_pos_list = tor_pos(nphi=4*n_plane) !###
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call apply_four_filter(result, simple_filter(m=0,n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    write(comment,'(a,i6.6)') 'time step #', index_start
    
    call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                     &
      filename=trim(filename), append=(.not.first_step), blanks=.true., comment=trim(comment))
    
  end subroutine average
  
  
  
  
  
  !> Output equilibrium information.
  subroutine equil_params(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, i_file
    character(len=1024) :: filename, status, access
    real*8 :: time
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    
    write(filename,'(4a)') DIR, 'equil_params',                                                    &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133
    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    if ( first_step ) then
      write(i_file,'(a)') '# time                R_axis              Z_axis              '//       &
        'Psi_axis            R_xpoint(1)         Z_xpoint(1)         Psi_xpoint(1)       '//       &
        'R_xpoint(2)         Z_xpoint(2)         Psi_xpoint(2)       R_lim               '//       &
        'Z_lim               Psi_lim             Psi_bnd'
    end if
    
    ! (not elegant, admittedly... but guarantees consistent time normalization:)
    call eval_expr(eq, units, exprs('t',1),                                                        &
      pol_pos(node_list,element_list,eq,R=eq%R_axis,Z=eq%Z_axis), tor_pos(phi=0.d0), result, ierr)
    time = result(1,1,1,1)
    
    write(i_file,'(es20.13,33f20.16)') time, eq%R_axis, eq%Z_axis, eq%Psi_axis, eq%R_xpoint(1),    &
      eq%Z_xpoint(1), eq%Psi_xpoint(1), eq%R_xpoint(2), eq%Z_xpoint(2), eq%Psi_xpoint(2), eq%R_lim,&
      eq%Z_lim, eq%Psi_lim, eq%Psi_bnd
    
    close(i_file)
    
  end subroutine equil_params
  
  
  
  
  
  !> Output the q-profile as a function of Psi_N
  recursive subroutine qprofile(command, first_step, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag

    
    ! --- Local variables
    integer                  :: k, k2, npts
    real*8, allocatable      :: q(:), rad(:)
    type (type_surface_list) :: surface_list
    character(len=1024)      :: filename, comment
    type(t_expr_list)        :: tmp_expr_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename,'(4a)') DIR, 'qprofile', trim(step_range_string(loop_min_step,loop_max_step)),  &
      '.dat'
    
    ! --- Find flux surfaces and determine q-profile
    surface_list%n_psi = npts
    allocate( surface_list%psi_values(npts), q(npts), rad(npts) )
    do k = 1, npts
      surface_list%psi_values(k) = eq%psi_axis + (eq%psi_bnd - eq%psi_axis) * real(k-1)/real(npts-1)
    end do
    call find_flux_surfaces(xpoint, xcase, node_list, element_list, surface_list)
    call determine_q_profile(node_list, element_list, surface_list, eq%psi_axis, eq%psi_xpoint,    &
      eq%Z_xpoint, q, rad)
    
    ! --- Write out q-profile versus Psi_n
    tmp_expr_list%n_expr = 0
    tmp_expr_list%expr(1)%name = 'Psi_n'
    tmp_expr_list%expr(2)%name = 'q'
    write(comment,'(a,i6.6)') 'time step #', index_start
    allocate(res1d(npts-2,2))
    do k2 = 1, npts-2
      k = k2 + 1 ! to avoid first and last point of q-profile which often is bad
      res1d(k2,:) = (/ get_psi_n(eq, surface_list%psi_values(k)), q(k) /)
    end do
    
    call write_ascii_1d(ierr, eq, tmp_expr_list, res1d, FORM_TABLE, header=.true.,                 &
      filename=trim(filename), append=(.not.first_step), blanks=.true., comment=trim(comment))
    
    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    if ( allocated(q)                          ) deallocate(q)
    if ( allocated(rad)                        ) deallocate(rad)
    
  end subroutine qprofile
  
  
  
  
  
  !> Output the flux surfaces.
  recursive subroutine fluxsurfaces(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    
    ! --- Local variables
    integer                  :: i, j, i_elm, npts, ip, nplot, i_file
    type (type_surface_list) :: surface_list
    character(len=1024)      :: filename, comment
    type(t_expr_list)        :: tmp_expr_list
    real*8                   :: psi_min, psi_max, psi_min2, psi_max2, ss1, dss1, ss2, dss2, tt1,   &
      dtt1, tt2, dtt2, u, si, dsi, ti, dti, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss,&
      Z_tt
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename,'(4a)') DIR, 'fluxsurfaces', trim(step_range_string(index_start,index_start)),  &
      '.dat'
    
    ! --- Find minimum and maximum psi values
    psi_min=+1.d99
    psi_max=-1.d99
    do i_elm = 1, element_list%n_elements
      call psi_minmax(node_list, element_list, i_elm, psi_min2, psi_max2)
      psi_min = min(psi_min,psi_min2)
      psi_max = max(psi_max,psi_max2)
    end do
    psi_min = psi_min + 0.001*(psi_max-psi_min)
    psi_max = psi_max - 0.001*(psi_max-psi_min)
    
    ! --- Find flux surfaces
    surface_list%n_psi = npts
    allocate( surface_list%psi_values(npts) )
    do i = 1, npts
      surface_list%psi_values(i) = psi_min + (psi_max-psi_min) * real(i-1)/real(npts-1)
    end do
    call find_flux_surfaces(xpoint, xcase, node_list, element_list, surface_list)
    
    ! --- Write out flux surfaces
    nplot  = 5
    i_file = 111
    call open_ascii_file(ierr, i_file, filename, .false.)
    do i = 1, npts
      
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
          write(i_file,'(2ES16.7)') R, Z
        end do
        
        write(i_file,*)
        write(i_file,*)
      
      end do
      
    end do
    
    close(i_file)
    
    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    
  end subroutine fluxsurfaces
  
  
  
  
  
  !> Output the separatrix.
  recursive subroutine separatrix(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer                  :: i, j, i_elm, npts, ip, nplot, i_file
    type (type_surface_list) :: surface_list
    character(len=1024)      :: filename, comment
    type(t_expr_list)        :: tmp_expr_list
    real*8                   :: psi_min, psi_max, psi_min2, psi_max2, ss1, dss1, ss2, dss2, tt1,   &
      dtt1, tt2, dtt2, u, si, dsi, ti, dti, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss,&
      Z_tt
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    write(filename,'(4a)') DIR, 'separatrix', trim(step_range_string(index_start,index_start)),    &
      '.dat'
    
    ! --- Find flux surfaces
    npts = 1
    surface_list%n_psi = 1
    allocate( surface_list%psi_values(1) )
    surface_list%psi_values(1) = eq%psi_bnd
    call find_flux_surfaces(xpoint, xcase, node_list, element_list, surface_list)
    
    ! --- Write out flux surfaces
    nplot  = 5
    i_file = 111
    call open_ascii_file(ierr, i_file, filename, .false.)
    i = 1
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
        call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt,      &
          Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
          
        ! --- Write out the (R,Z)-coordinates
        write(i_file,'(2ES16.7)') R, Z
      end do
      
      write(i_file,*)
      write(i_file,*)
    
    end do
    
    close(i_file)
    
    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    
  end subroutine separatrix
  
  
  
  
  
  !> Perform a 2D Fourier analysis (in straight field line coordinates).
  subroutine four2d(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts
    character(len=1024) :: filename_start
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename_start,'(3a)') DIR, 'exprs_four2d', trim(step_range_string(index_now,index_now))
    
    call fourier_analysis(node_list, element_list, eq, units, expr_list, cp, npts, ierr,           &
      filename_start, OUTP_ABS_VALUE)
    
  end subroutine four2d
  
  
  
  
  
  !> Output magnetic field for the Gourdon code.
  subroutine gourdon(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts, n_R, n_Z, n_phi, i, j, k
    character(len=1024) :: filename
    real*8 :: R_min, R_max, Z_min, Z_max, R_max2, Z_max2, phi_max2, fact_phi, fact_btor, fact_bpol,&
      tmp
    real*8, allocatable :: field(:,:,:,:)
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    type(t_expr_list),    save :: tmp_expr_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,10); if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    ! --- Preparation
    R_min = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    R_max = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    n_R   = to_int  (command%args(3), ierr); if ( ierr /= 0 ) return
    Z_min = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    Z_max = to_float(command%args(5), ierr); if ( ierr /= 0 ) return
    n_Z   = to_int  (command%args(6), ierr); if ( ierr /= 0 ) return
    n_phi = to_int  (command%args(7), ierr); if ( ierr /= 0 ) return
    !   --- The following three parameters are correction factors which allow to
    !     transform the field between different coordinate systems.
    fact_phi  = to_float(command%args(8),  ierr); if ( ierr /= 0 ) return
    fact_btor = to_float(command%args(9),  ierr); if ( ierr /= 0 ) return
    fact_bpol = to_float(command%args(10), ierr); if ( ierr /= 0 ) return
    
    write(filename,'(3a)') DIR, 'gourdon', trim(step_range_string(index_now,index_now))
    
    ! --- Take into account that the last points are not included in Gourdon format!
    R_max2   = R_min + real(n_R-1)/real(n_R) * (R_max-R_min)
    Z_max2   = Z_min + real(n_Z-1)/real(n_Z) * (Z_max-Z_min)
    phi_max2 = real(n_phi-1)/real(n_phi) * 2.d0 * pi         * fact_phi
    
    ! --- Make sure that positions outside the JOREK domain get bfield=0.
    tmp = expr_outside_value
    expr_outside_value = 0.d0
    
    ! --- Calculate field components.
    if ( first_step ) then ! (Positions remain unchanged for all time steps, compute only once)
      call create_pol_pos(pol_pos_list, ierr, node_list, element_list, eq, Rmin=R_min, Rmax=R_max2,&
        nR=n_R-1, Zmin=Z_min, Zmax=Z_max-1, nZ=n_Z-1)
      tor_pos_list  = tor_pos(phistart=0.d0, phiend=phi_max2, nphi=n_phi-1)
      tmp_expr_list = exprs((/'B_tor', 'B_R  ', 'B_Z  '/), 3)
    end if
    call eval_expr(eq, JOREK_UNITS, tmp_expr_list, pol_pos_list, tor_pos_list, result, ierr)
    if ( fact_btor /= 1.d0 ) result(:,:,:,1  ) = result(:,:,:,1  ) * fact_btor
    if ( fact_bpol /= 1.d0 ) result(:,:,:,2:3) = result(:,:,:,2:3) * fact_bpol
    allocate(field(3,n_R-1,n_Z-1,n_phi-1))
    !    --- Array transform for file output...
    do k = 1, n_Z-1
      do j = 1, n_R-1
        do i = 1, 3
          field(i,j,k,:) = result(:,j,k,i)
        end do
      end do
    end do
    
    ! --- Write out the data.
    open(122, file=filename, status='replace', action='write', form='unformatted', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) 'ERROR in routine gourdon: Creating file "', trim(filename), '" failed.'
      return
    end if
    write(122) real(n_phi), real(n_R), real(n_Z), real(3), real(1), real(n_phi)
    write(122) (R_max+R_min)/2.d0, (Z_max+Z_min)/2.d0, (R_max-R_min)/2.d0, (Z_max-Z_min)/2.d0
    do k = 0, n_phi - 1
      write(122) field(:,:,:,k)
    end do
    close(122)
    
    ! --- Clean up.
    expr_outside_value = tmp ! restore previous value
    if ( allocated(field ) ) deallocate(field)
    
  end subroutine gourdon
  
  
  
  
  
  !> Select JOREK normalized units.
  subroutine select_jorek_units(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0); if ( ierr /= 0 ) return
    
    call set_setting('units', '0', ierr)
    
  end subroutine select_jorek_units
  
  
  
  
  
  !> Select SI units.
  subroutine select_si_units(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0); if ( ierr /= 0 ) return
    
    call set_setting('units', '1', ierr)
    
  end subroutine select_si_units
  
  
  
  
  
end module exec_commands

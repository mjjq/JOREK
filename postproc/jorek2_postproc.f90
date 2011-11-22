!> Program for post-processing data written out as JOREK restart files (jorekXXXXX.rst).
!!
!! * Interactive use: Call './jorek2_postproc' in your run folder. You will get some
!!   'getting started' information. Enter 'help' to see a list of available commands
!!   or 'help <command>' for command-specific usage information.
!! * Script: Run './jorek2_postproc < script'
program jorek2_postproc
  
  use nodes_elements, only: node_list, element_list
  use phys_module
  use parse_commands, only: read_command, type_command
  use exec_commands,  only: exec_command, specific_help
  use settings,       only: set_setting
  
  implicit none
  
  integer, parameter :: file_handle = 17
  type(type_command) :: command
  integer            :: error
  
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Preset some parameters
  call set_setting('linepoints',      '200',   error)
  call set_setting('tor_points',      '50',    error)
  call set_setting('surfaces',        '100',   error)
  call set_setting('verbose',         'false', error)
  call set_setting('debug',           'false', error)
  
  ! --- Print getting started information
  call specific_help('getting_started')
  
  do ! (main loop: Read, parse, and execute one command after the other)
    
    ! --- Read and parse a command line
    call read_command(command, error)
    if ( error /= 0 ) then
      write(*,*) 'ERROR in read_command'
      cycle
    else if ( command%n_options == 0 ) then
      cycle
    end if
    
    ! --- Exit upon user request
    if ( (command%option(1) == 'exit') .or. (command%option(1) == 'quit') ) stop
    
    ! --- Execute the command
    call exec_command(command, .true., file_handle, error)
    if ( error /= 0 ) then
      write(*,*) 'ERROR in exec_command'
      cycle
    end if
    
  end do ! (main loop)
  
end program jorek2_postproc

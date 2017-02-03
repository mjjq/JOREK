!> Module for reading and parsing a command (used by jorek2_postproc)
module parse_commands

  use convert_character,    only: lower_case
  use settings,             only: set_setting, get_setting
  
  implicit none
  
  integer, parameter :: MAX_ARGS    = 100 !< Maximum number of arguments for a command
  integer, parameter :: FILE_HANDLE = 16  !< File handle
  
  integer, parameter :: STDIN_MODE  = 1   !< Commands are read from STDIN in this mode
  integer, parameter :: SOURCE_MODE = 2   !< Commands are read from a file in this mode
  integer            :: mode = STDIN_MODE !< Current mode (STDIN_MODE or SOURCE_MODE)
  
  !> Type declaration for a command entered by the user
  !!   arg(0) is the command itself
  !!   arg(1:n_args) are the command parameters
  type :: type_command
    integer            :: n_args           ! Number of arguments
    character(len=256) :: args(0:MAX_ARGS) ! Actual arguments; args(0) is the command itself
  end type type_command
  
  
  
  private
  public type_command, print_command, read_command
  
  
  
  save
  
  
  
  contains
  
  
  
  !> Read a command
  recursive subroutine read_command(command, error)
    
    ! --- Routine parameters
    type(type_command), intent(inout) :: command !< The command read from STDIN is returned
    integer,            intent(out)   :: error   !< Error flag
    
    ! --- Local variables
    character(len=1024) :: line
    integer             :: ierr
    
    error = 0
    
    ! --- Read next command line
    if ( mode == STDIN_MODE ) then ! ... from STDIN
      write(*,'(a)',advance='no') '> '
      read(*,'(a)',iostat=ierr) line
      if ( ierr /= 0 ) line = 'exit' ! End of input (Ctrl-D) is equivalent to 'exit' command
    else ! ... from a file
      read(FILE_HANDLE,'(a)',iostat=ierr) line
      if ( ierr /= 0 ) then
        if ( get_setting('debug',error) == 'true' ) write(*,*) 'End of source file -- switch to interactive mode'
        close(FILE_HANDLE)
        mode = STDIN_MODE
        call read_command(command, error)
        return
      end if
    end if
    
    ! --- Parse the line
    call parse_command(line, command, error)
    if ( error /= 0 ) then
      write(*,*) 'ERROR in read_command executing parse_command.'
      return
    end if
    if ( get_setting('debug',error) == 'true' ) then
      write(*,'(a,a,a)') 'Entered: "', trim(line), '"'
      call print_command(command)
      write(*,*)
    end if
    
    ! --- If a script file is sourced, open it and switch to SOURCE_MODE
    if ( command%args(0) == 'source' ) then
      if ( command%n_args /= 1 ) then
        write(*,*) 'ERROR: "source" called with wrong number of arguments.'
        call read_command(command, error)
        return
      end if
      open(FILE_HANDLE, file=trim(command%args(1)), action='read', status='old', iostat=ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR: The script file "', trim(command%args(1)), '" could not be found.'
        call read_command(command, error)
        return
      end if
      mode = SOURCE_MODE
      call read_command(command, error)
      return
    end if
    
  end subroutine read_command
  
  
  
  !> Parse a command
  subroutine parse_command(line, command, error)
    
    ! --- Routine parameters
    character(len=*),   intent(in)    :: line    !< Character string containing the command
    type(type_command), intent(inout) :: command !< Parsed command
    integer,            intent(out)   :: error   !< Error flag
    
    ! --- Local variables
    character(len=LEN(line))        :: todo_string ! Copy of line
    character(len=1)                :: quote       ! Which quote (' or ")?
    integer                         :: iend        ! End position of argument in todo_string

    error = 0
    todo_string = trim(adjustl(line))
    command%n_args  = -1
    command%args(:) = ''
    iend = 0
    
    ! --- Extract one argument after the other from the command line
    do
        
      if (LEN_TRIM(todo_string) == 0) exit
      
      ! --- Argument in quotes
      if ( (todo_string(1:1) == '"') .or. (todo_string(1:1) == '''') ) then
        
        quote = todo_string(1:1)
        todo_string =  adjustl(todo_string(2:LEN(todo_string)))
        iend = INDEX(todo_string, quote) - 1
        if (iend == -1) then 
          write(*,*) 'ERROR: No matching quote (',quote,') found.'
          error = 1
          return
        else
          command%n_args = command%n_args + 1
          command%args(command%n_args) = todo_string(1:iend)
          todo_string =  adjustl(todo_string(iend+2:LEN(todo_string)))
        end if
      ! --- Comment
      else if (todo_string(1:1) == '#') then
        exit
      ! --- Normal argument
      else
        iend = INDEX(todo_string, ' ') - 1
        command%n_args = command%n_args + 1
        command%args(command%n_args) = todo_string(1:iend)
        todo_string =  adjustl(todo_string(iend+2:LEN(todo_string)))
      end if
    end do
    
  end subroutine parse_command
  
  
  
  !> Print a command (for debugging purposes)
  subroutine print_command(command)
    
    ! --- Routine parameters
    type(type_command), intent(in) :: command !< Command to be printed
    
    integer :: i
    
    if ( command%n_args >= 0 ) then
      write(*,'(1x,a,i3,a)',advance='no') 'Command: "'
      do i = 0, command%n_args
        if ( i /= 0 ) write(*,'(a)',advance='no') ' '
        write(*,'(a)',advance='no') trim(command%args(i))
      end do
      write(*,*) '"'
    else
      write(*,'(a)') '--EMPTY COMMAND--'
    end if
    
  end subroutine print_command
  
  
  
  !> Routine to run parse_command with some strange test cases.
  subroutine test_parse_command()
    type(type_command) :: command
    integer :: i, error
    character(len=60), parameter :: test_string(15) = (/                 &
      '                                                            ',    &
      ' ''                                                          ',   &
      ' ''''                                                         ',  &
      ' "                                                          ',    &
      ' ""                                                         ',    &
      'abc de f ghijklmnopqrstuvwxyz ''abra cadabra'' " blub blub "  ',  &
      'ifort -c -g -check all -debug all -warn all -ftrapuv a.f90  ',    &
      'abra                                                 cadabra',    &
      '''abra                                               cadabra''',  &
      '"abra                                               cadabra"',    &
      '#                                                           ',    &
      ' ## ## abra cadabra                                         ',    &
      ' #''"                                                        ',    &
      '"abra" "cadabra"                                            ',    &
      'line temperature 1.6 0.0 0.0 1.2 0.0 0.0                    ' /)
    
    call set_setting('debug',   'true', error)
    call set_setting('verbose', 'true', error)
    do i = 1, size(test_string,1)
      write(*,*)
      write(*,*)
      write(*,'(3a)') 'line=<', test_string(i), '>'
      call parse_command(test_string(i), command, error)
      if ( error /= 0 ) write(*,'(a)') '--an error occured--'
      call print_command(command)
    end do
    
  end subroutine test_parse_command
  
end module parse_commands

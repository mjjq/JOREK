!> Module containing the usage help texts
!! (used by jorek2_postproc)
module postproc_help
  
  
  
  implicit none
  
  
  
  contains
  
  
  
  !> Print general usage information
  subroutine general_help()
    
    write(*,*) ''
    write(*,*) '============================================================================'
    write(*,*) 'Help topics:'
    write(*,*) '  getting_started   Very brief introduction by an example'
    write(*,*) '  average           Poloidally and toroidally averaged profiles'
    write(*,*) '  exit / quit       Exit the program'
    write(*,*) '  equil_params      Output equilibrim parameters (axis, X-point, limiter)'
    write(*,*) '  expressions       List or select physical expressions'
    write(*,*) '  for               Loop over one or several timesteps'
    write(*,*) '  mark_coords       Mark expressions as coordinates'
    write(*,*) '  midplane          Toroidally averaged expressions on the midplane'
    write(*,*) '  namelist          Load a namelist input file'
    write(*,*) '  params            Print JOREK parameters'
    write(*,*) '  point             Expressions at a certain (R,Z,phi) position'
    write(*,*) '  pol_line          Expressions along line from (R0,Z0,phi) to (R1,Z1,phi)'
    write(*,*) '  qprofile          Export q-profile versus Psi_N'
    write(*,*) '  set               Change a setting or list all postproc settings'
    write(*,*) '  timesteps         List available restart files'
    write(*,*) '  tor_line          Expressions along line from (R,Z,phi0) to (R,Z,phi1)'
    write(*,*) '----------------------------------------------------------------------------'
    write(*,*) '  Enter "help <command>" for details'
    write(*,*) '============================================================================'
    write(*,*) ''
    
  end subroutine general_help
  
  
  
  !> Print specific help for a certain command
  subroutine specific_help(topic)
    
    character(len=*), intent(in) :: topic !< Print specific help for this help topic
    
    write(*,*) ''
    write(*,*) '========================================================================='
    write(*,*) 'Help topic: ', trim(topic)
    write(*,*) '-------------------------------------------------------------------------'
    select case ( trim(topic) )
      case ( 'getting_started' )
        write(*,*) 'Welcome to the interactive postprocessing tool JOREK2_POSTPROC!'
        write(*,*) ''
        write(*,*) 'An example to get you started:'
        write(*,*) ''
        write(*,*) '  namelist input               # Read namelist input file'
        write(*,*) '  for step 100 to 200 do       # Do the following for severl time steps:'
        write(*,*) ''
        write(*,*) '    expressions Psi_N T rho      # Select physical expressions'
        write(*,*) '    mark_coords 1                # Mark first expression as coordinate'
        write(*,*) '    average                      # Toroidal and poloidal average'
        write(*,*) ''
        write(*,*) '    expressions R rho T Psi Er   # Select expressions'
        write(*,*) '    mark_coords 1                # Mark first expression as coordinate'
        write(*,*) '    midplane                     # Midplane profiles'
        write(*,*) '  done'
        write(*,*) ''
        write(*,*) 'All output is written into subfolder postproc/'
        write(*,*) ''
        write(*,*) 'Type "help" for general usage information or, e.g., "help average" for'
        write(*,*) 'command-specific help.'
        write(*,*) ''
        write(*,*) 'A documentation is available at http://jorek.eu/wiki.'
      case ( 'average' )
        write(*,*) 'Usage:'
        write(*,*) '  average'
        write(*,*) ''
        write(*,*) 'Determine poloidally and toroidally averaged profiles.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set surfaces 200'
        write(*,*) '  expressions Psi_N rho T Er'
        write(*,*) '  mark_coords 1'
        write(*,*) '  for step 1 to 200'
        write(*,*) '    average'
        write(*,*) '  done'
      case ( 'equil_params' )
        write(*,*) 'Usage:'
        write(*,*) '  equil_params'
        write(*,*) ''
        write(*,*) 'Output equilibrim parameters (axis, X-point, limiter).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 1 to 200'
        write(*,*) '    equil_params'
        write(*,*) '  done'
      case ( 'expressions' )
        write(*,*) 'Usage:'
        write(*,*) '  expressions [<expr1> <expr2> ...]'
        write(*,*) ''
        write(*,*) 'Select physical expressions. Required by many commands. If called without'
        write(*,*) 'parameters, a table of available expressions is printed.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  expressions'
        write(*,*) '  expressions Psi_N rho T Er'
      case ( 'for' ) 
        write(*,*) 'Usage:'
        write(*,*) '  for step <from> [to <to>] do'
        write(*,*) '    ...'
        write(*,*) '  done'
        write(*,*) ''
        write(*,*) 'Execute one or several commands for one or several time-steps.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 1 do'
        write(*,*) '    ...'
        write(*,*) '  done'
        write(*,*) ''
        write(*,*) '  for step 1 to 200 do'
        write(*,*) '    ...'
        write(*,*) '  done'
      case ( 'mark_coords' )
        write(*,*) 'Usage:'
        write(*,*) '  mark_coords <n>'
        write(*,*) ''
        write(*,*) 'Mark the first <n> expressions as coordinates. Required by some commands.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  expressions Psi_N rho T Er'
        write(*,*) '  mark_coords 1'
      case ( 'midplane' )
        write(*,*) 'Usage:'
        write(*,*) '  midplane'
        write(*,*) ''
        write(*,*) 'Output toroidally averaged midplane profiles.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 200'
        write(*,*) '  for step 1 to 100 do'
        write(*,*) '    expressions R rho T Psi Er'
        write(*,*) '    mark_coords 1                    # mark first expression as coordinate'
        write(*,*) '    midplane'
        write(*,*) '  done'
      case ( 'namelist' )
        write(*,*) 'Usage:'
        write(*,*) '  namelist <namelist input file>'
        write(*,*) ''
        write(*,*) 'Load namelist input file. Required by most commands.'
        write(*,*) ''
        write(*,*) 'Example:'
        write(*,*) '  namelist input1'
      case ( 'params' )
        write(*,*) 'Usage:'
        write(*,*) '  params'
        write(*,*) ''
        write(*,*) 'Output JOREK parameters (hardcoded and namelist input).'
        write(*,*) ''
      case ( 'point' )
        write(*,*) 'Usage:'
        write(*,*) '  point <R> <Z> <phi>'
        write(*,*) ''
        write(*,*) 'Determine the value of a variable at position (R,Z,phi).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 1 to 100 do'
        write(*,*) '    expressions rho T Psi Er'
        write(*,*) '    point 1.6 0.0 0.0'
        write(*,*) '  done'
      case ( 'pol_line' )
        write(*,*) 'Usage:'
        write(*,*) '  pol_line <R0> <Z0> <R1> <Z1> <phi>'
        write(*,*) ''
        write(*,*) 'Expressions along line from (R0,Z0,phi) to (R1,Z1,phi).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 500'
        write(*,*) '  expressions R rho T E_r'
        write(*,*) '  mark_coords 1'
        write(*,*) '  for step 0 to 100 do'
        write(*,*) '    pol_line 1.6 0.5 1.7 0.5 0.'
        write(*,*) '  done'
      case ( 'set' )
        write(*,*) 'Usage:'
        write(*,*) '  set [<name> <value>]'
        write(*,*) ''
        write(*,*) 'Sets a setting to a certain value. Lists all settings if called'
        write(*,*) 'without parameters.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 500     # used by pol_line, midplane, ...'
        write(*,*) '  set surfaces 200       # used by fluxsurfaces, average, ...'
        write(*,*) '  set units 0            # 0: JOREK units; 1: SI units'
        write(*,*) '  set'
      case ( 'timesteps' )
        write(*,*) 'Usage:'
        write(*,*) '  timesteps'
        write(*,*) ''
        write(*,*) 'Lists all available time steps (restart files).'
      case ( 'tor_line' )
        write(*,*) 'Usage:'
        write(*,*) '  tor_line <R> <Z> <phi0> <phi1>'
        write(*,*) ''
        write(*,*) 'Expressions along line from (R,Z,phi0) to (R,Z,phi1).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 500'
        write(*,*) '  expressions phi rho T E_r'
        write(*,*) '  mark_coords 1'
        write(*,*) '  for step 0 to 100 do'
        write(*,*) '    tor_line 1.6 0.1 0. 3.141592'
        write(*,*) '  done'
      case default
        write(*,*) '[This help topic does not exist, sorry]'
    end select
    write(*,*) '========================================================================='
    write(*,*) ''
    
  end subroutine specific_help
  
  
  
end module postproc_help

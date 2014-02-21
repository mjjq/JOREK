!> Module containing the usage help texts
!! (used by jorek2_postproc)
module postproc_help
  
  
  
  implicit none
  
  
  
  contains
  
  
  
  !> Print general usage information
  subroutine general_help()
    
    write(*,*) ''
    write(*,*) '============================================================================'
    write(*,*) 'Available commands:'
    write(*,*) '  average           Poloidally and toroidally averaged profiles'
    write(*,*) '  exit / quit       Exit the program'
    write(*,*) '  expressions       List or select expressions'
    write(*,*) '  for ...           Loop over one or several timesteps'
!    write(*,*) '  fluxsurfaces      Write out flux surfaces'
!    write(*,*) '  global_parameters Output global parameters (beta_p, plasma current, ...)'
!    write(*,*) '  gourdon           Export the magnetic field for Gourdon'
!    write(*,*) '  heatfluxpattern   Determine the heat flux pattern on a target plate'
!    write(*,*) '  line              Variable-value along a line in (R,Z,phi)'
    write(*,*) '  mark_coords       Mark expressions as coordinates'
    write(*,*) '  midplane          Toroidally averaged expressions on the midplane'
    write(*,*) '  namelist          Load a namelist input file'
    write(*,*) '  params            Print JOREK parameters'
    write(*,*) '  point             Expressions at a certain (R,Z,phi) position'
!    write(*,*) '  qprofile          Output the q-profile'
    write(*,*) '  set               Change a setting or list all settings'
    write(*,*) '  timesteps         List available restart files'
!    write(*,*) '  volume            Calculate the plasma volume'
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
    write(*,*) '###help not yet updated###'
    return
    select case ( trim(topic) )
      case ( 'getting_started' )
        write(*,*) 'Welcome to the interactive postprocessing tool JOREK2_POSTPROC!'
        write(*,*) ''
        write(*,*) 'To get started, you need to load a JOREK namelist input file, e.g.:'
        write(*,*) '  namelist input'
        write(*,*) ''
        write(*,*) 'Then, execute commands for one or several timesteps, e.g.:'
        write(*,*) '  for step 100 to 200 do'
        write(*,*) '    average temperature'
        write(*,*) '    fluxsurfaces'
        write(*,*) '    point density 1.6 0 0'
        write(*,*) '  done'
        write(*,*) ''
        write(*,*) 'Type "help" for general usage information or, e.g., "help average" for'
        write(*,*) 'command-specific help.'
      case ( 'average' )
        write(*,*) 'Usage:'
        write(*,*) '  average <variable> [<variable> ...]'
        write(*,*) ''
        write(*,*) 'Determine the flux-surface average of a variable.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set surfaces 200'
        write(*,*) '  for step 1 to 200'
        write(*,*) '    average temperature density flux'
        write(*,*) '  done'
      case ( 'axis' )
        write(*,*) 'Usage:'
        write(*,*) '  axis'
        write(*,*) ''
        write(*,*) 'Locate the magnetic axis over multiple timesteps'
        write(*,*) ''
        write(*,*) 'Examples:'        
        write(*,*) '  for step 100 to 200 do'
        write(*,*) '    axis'
        write(*,*) '  done'
        write(*,*) ''
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
        write(*,*) '    fluxsurfaces'
        write(*,*) '  done'
        write(*,*) ''
        write(*,*) '  for step 1 to 200 do'
        write(*,*) '    average temperature'
        write(*,*) '    line flux 1.6 0.0 0.0 1.7 0.0 0.0'
        write(*,*) '    point density 1.6 0.0 0.0'
        write(*,*) '  done'
      case ( 'fluxsurfaces' ) 
        write(*,*) 'Usage:'
        write(*,*) '  fluxsurfaces'
        write(*,*) ''
        write(*,*) 'Determine the flux surfaces.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set surfaces 200'
        write(*,*) '  for step 550 do'
        write(*,*) '    fluxsurfaces'
        write(*,*) '  done'
      case ( 'gourdon' ) 
        write(*,*) 'Usage:'
        write(*,*) '  gourdon <R_min> <R_max> <n_R> <Z_min> <Z_max> <n_Z> <n_phi> [<derivs>]'
        write(*,*) ''
        write(*,*) 'Export the magnetic field for the Gourdon code.'
        write(*,*) ''
        write(*,*) 'Details:'
        write(*,*) '  <R_min>, <R_max>:      R-range to export the field in'
        write(*,*) '  <Z_min>, <Z_max>:      Z-range -"-'
        write(*,*) '  <n_R>, <n_Z>, <n_phi>: Number of points in R-, Z-, phi-directions'
        write(*,*) '  <derivs>:              Include field derivatives? Can be .true. or'
        write(*,*) '                         .false. (default).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 550 do'
        write(*,*) '    gourdon 1.5 1.7 200 -0.1 0.1 200 200 .true.'
        write(*,*) '  done'
      case ( 'global_parameters' ) 
        write(*,*) 'Usage:'
        write(*,*) '  global_parameters'
        write(*,*) ''
        write(*,*) 'Output global parameters (beta_p, plasma current, ...).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 550 do'
        write(*,*) '    global_parameters'
        write(*,*) '  done'
      case ( 'heatfluxpattern' ) 
        write(*,*) 'Usage:'
        write(*,*) '  heatfluxpattern <target-plate-file> <n_0> <m_ion> [all|diffusive|convective]'
        write(*,*) ''
        write(*,*) 'Determine the heat flux pattern at an axis-symmetric target plate.'
        write(*,*) ''
        write(*,*) 'Details:'
        write(*,*) '  <target-plate-file>: Ascii file containing line segments in (R,Z)'
        write(*,*) '  <n_0>:               Core ion density [m^-3]'
        write(*,*) '  <m_ion>:             Mass of an ion [kg]; m_D=3.43e-27'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 550 do'
        write(*,*) '    heatfluxpattern target-plate.dat 8.e19 3.343e-27'
        write(*,*) '  done'
      case ( 'lfs_hfs' )
        write(*,*) 'Usage:'
        write(*,*) '  lfs_hfs <variable>'
        write(*,*) ''
        write(*,*) 'Output toroidally averaged low and high field side profiles of a variable.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 500'
        write(*,*) '  for step 550 do'
        write(*,*) '    lfs_hfs temperature'
        write(*,*) '  done'
      case ( 'line' )
        write(*,*) 'Usage:'
        write(*,*) '  line <variable1> <variable2> <R0> <Z0> <phi0> <R1> <Z1> <phi1>'
        write(*,*) ''
        write(*,*) 'Determine the value of variable2 versus variable1 along a line in (R,Z,phi).'
        write(*,*) ''
        write(*,*) 'Valid values for variable1 and variable 2 are all JOREK variables and'
        write(*,*) '"R", "Z", "phi", "normalized_length", and "Psi_N".'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 500'
        write(*,*) '  for step 550 do'
        write(*,*) '    line flux normalized_length 1.6 0.0 0.0 1.7 0.0 0.0'
        write(*,*) '    line density flux 1.7 0.1 0.2 1.8 0.2 0.3'
        write(*,*) '  done'
      case ( 'namelist' )
        write(*,*) 'Usage:'
        write(*,*) '  namelist <namelist input file>'
        write(*,*) ''
        write(*,*) 'Load namelist input file'
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
        write(*,*) '  point <variable> <R> <Z> <phi>'
        write(*,*) ''
        write(*,*) 'Determine the value of a variable at position (R,Z,phi).'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  point flux 1.6 0.0 0.0'
        write(*,*) '  point 3    1.7 0.1 0.2'
      case ( 'qprofile' )
        write(*,*) 'Usage:'
        write(*,*) '  qprofile'
        write(*,*) ''
        write(*,*) 'Output the q-profile as a function of Psi_N.'
        write(*,*) ''
      case ( 'set' )
        write(*,*) 'Usage:'
        write(*,*) '  set [<name> <value>]'
        write(*,*) ''
        write(*,*) 'Sets a setting to a certain value. Lists all settings if called'
        write(*,*) 'without parameters.'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  set linepoints 500'
        write(*,*) '  set surfaces 200'
        write(*,*) '  set'
      case ( 'volume' )
        write(*,*) 'Usage:'
        write(*,*) '  volume'
        write(*,*) ''
        write(*,*) 'Calculate the plasma volume'
        write(*,*) ''
        write(*,*) 'Examples:'
        write(*,*) '  for step 500 do'
        write(*,*) '    volume'
        write(*,*) '  done'
      
      case default
        write(*,*) '[This help topic does not exist, sorry]'
    end select
    write(*,*) '========================================================================='
    write(*,*) ''
    
  end subroutine specific_help
  
  
  
end module postproc_help

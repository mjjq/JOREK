!> Program for post-processing jorek data.
!!
!! Documentation at http://jorek.eu/wiki
program jorek2_postproc
  
  use nodes_elements, only: node_list, element_list
  use phys_module
  use mod_new_diag
  use parse_commands, only: read_command, type_command
  use exec_commands,  only: exec_command, specific_help
  use settings,       only: set_setting
  use basis_at_gaussian, only: initialise_basis
  use mpi_mod
  
  implicit none
  
  type(type_command) :: command
  integer            :: ierr
  integer            :: i_file, i, units, my_id, n_cpu
  integer            :: required, provided, StatInfo
  integer*4          :: rank, comm_size 

  ! --- Initialise MPI / threaded MPI
#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
#ifdef STAN_FLAG
required = 0
#endif
  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
  
  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size

  ! --- jorek2_postproc is not ready yet for MPI
  if (n_cpu /= 1) then
    write(*,*) "Please execute with mpirun -n 1 ./jorek2_postproc, multi MPI is not ready yet"
    stop
  endif
  
  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank
 
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Initialize the new_diag package.
  call init_new_diag(.false.)
  
  ! --- Preset namelist input parameters
  call preset_parameters()
  
  ! --- Preset some parameters
  call set_setting('units',           '0',     ierr) ! JOREK units ("set units 1" for SI-units)
  call set_setting('loop_units',      '0',     ierr) ! JOREK units ("set units 1" for SI-units)
  call set_setting('linepoints',      '200',   ierr)
  call set_setting('tor_points',      '50',    ierr)
  call set_setting('surfaces',        '100',   ierr)
  call set_setting('verbose',         'false', ierr)
  call set_setting('debug',           'false', ierr)
  call set_setting('nsmallsteps',     '3',     ierr) ! numerical parameter for straight field lines
  call set_setting('nmaxsteps',       '2500',  ierr)
  call set_setting('deltaphi',        '0.3',   ierr)
  call set_setting('rad_range_min',   '0.001', ierr)
  call set_setting('rad_range_max',   '0.999', ierr)
  call set_setting('nTht',            '32',    ierr)
  
  ! --- Print getting started information
  call specific_help('getting_started')
  
  do ! (main loop: Read, parse, and execute one command after the other)
    
    ! --- Read and parse a command line
    call read_command(command, ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR in read_command'
      cycle
    else if ( command%n_args < 0 ) then
      cycle
    end if
    
    ! --- Exit upon user request
    if ( (command%args(0) == 'exit') .or. (command%args(0) == 'quit') ) stop
    
    ! --- Execute the command
    call exec_command(command, .true., ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR in exec_command'
      cycle
    end if
    
  end do ! (main loop)
  
end program jorek2_postproc

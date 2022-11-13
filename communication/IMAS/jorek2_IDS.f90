!> Program to put JOREK data into IMAS IDSs
program jorek2_IDS

#ifdef USE_IMAS
  use ids_schemas 
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get, ids_put

  use mod_jorek2IMAS 
  use nodes_elements
  use data_structure
  use mod_import_restart
  use phys_module, only : rst_format
  use basis_at_gaussian, only: initialise_basis
  
  implicit none
  
  character(len=200):: user, database
  character(len=64) :: file_name, name_proj
  integer :: shot_number, run_number, i_begin, i_end, i_step
  integer :: ierr, idx, stat
  logical :: first_step, file_exists
  logical :: export_MHD, export_radiation, only_proj

  namelist /imas_params/ shot_number, run_number, user, database, i_begin, i_end, &
                         export_mhd, export_radiation, only_proj

  ! --- Necessary initialization ------------------
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Preset namelist input parameters
  call preset_parameters()

  ! --- Read input file
  call initialise_parameters(0, "__NO_FILENAME__")
  ! -----------------------------------------------
  
  ! --- Preset parameters for this program
  database    = 'test'                    ! Name of the database to export the results
  shot_number = 111112;   run_number=1;   
  i_begin     = 0                         ! Starting restart file index
  i_end       = 99999                     ! Ending restart file index
  export_MHD       = .true.
  export_radiation = .false.
  only_proj        = .false.              ! true if only projection*.h5 files are available
                                          ! a normal JOREK restart file must be copied into jorek_restart.h5

  call getenv('USER',user)
  
  ! --- Read parameters from namelist file 'imas.nml' if it exists
  open(42, file='imas.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    write(*,*) 'Reading parameters from imas.nml namelist.'
    read(42, imas_params)
    close(42)
  end if

  ! --- Try to open shot and number if it exists
  call imas_open_env( 'ids', shot_number,run_number,idx,user,database,'3',stat)! 3 is the database version  
  
  if (stat /= 0) then  ! --- Create a new shot if it doesn't exist
    write(*,*) '  Shot/run number did not exist, creating new one...'
    call imas_create_env('ids',shot_number,run_number, 0,0,idx,user,database,'3') 
  endif

  if (export_radiation)  allocate( aux_node_list )

  first_step = .true.

  ! --- Loop over
  do i_step = i_begin, i_end
 
    ! --- Cycle when required files don't exist 
    if (only_proj .and. export_radiation) then
      write(name_proj,'(a,i9.9,a)') 'projections000.', i_step, '.h5'  ! This formatting should be improved
      inquire (file=trim(name_proj), exist=file_exists)
      if (.not. file_exists) cycle
      file_name = 'jorek_restart' 
    else
      if (.not. restart_file_exists(i_step)) cycle
      write(file_name,'(a,i5.5)')   'jorek', i_step
      write(name_proj,'(a,i5.5,a)') 'projections', i_step, '.h5'
    endif

    ! --- Import restart file
    write(*,*)
    write(*,'(a,i9.9,a)') '#################### STEP ', i_step, ' ####################'
    write(*,*)
    call import_restart(node_list, element_list, file_name, rst_format, ierr)
    if (ierr /=0 ) then
       write(*,*) '  Could not read the restart file'
       stop
    endif

    ! --- Fill and export an MHD IDS
    if (export_mhd)  call fill_mhd_IDS(first_step, idx)  

    ! --- Fill and export a radiation IDS
    if (export_radiation) then
      call import_hdf5_restart_aux(aux_node_list, name_proj, rst_format, ierr)
      if (ierr /= 0) then
        write(*,*) ' Could not open projections file were radiation is stored'
        stop
      endif
      call fill_radiation_IDS(first_step, idx)  
    endif

    first_step = .false.

  enddo

  call imas_close(idx) 
 
#else

  write(*,*) 'Error: jorek2_IDS must be compiled with IMAS (USE_IMAS=1)'
  write(*,*) 'You will also have to load the IMAS module, in case you have not done it'
  write(*,*) '     module load IMAS                                                   '

#endif

end program jorek2_IDS

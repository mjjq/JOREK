!< This program reads a jorek_restart file and a "starwall_response.dat" file
!! and calculates the total wall forces on the STARWALL wall. It involves the 
!! computation of expensive volume integrals in the plasma, so run it
!! in paralllel with several MPI processes for large cases
program jorek2_wall_forces 

  use mod_vacuum_fields,   only: total_wall_forces 
  use data_structure
  use phys_module
  use mod_parameters
  use nodes_elements
  use mod_boundary,            only: boundary_from_grid
  use vacuum
  use vacuum_response,     only: get_vacuum_response, update_response, init_wall_currents, I_coils
  use vacuum_equilibrium,  only: import_external_fields
  use mod_import_restart
  use mod_element_rtree, only: populate_element_rtree
  use basis_at_gaussian, only: initialise_basis
  use tr_module
#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
  use matio_module, only: timestamp
#endif
  use mpi_mod

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  implicit none
 
  integer   :: my_id, my_id_n, my_id_master, ierr
  integer   :: i_rank(n_tor), n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer   :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  integer   :: required,provided,StatInfo
  integer*4 :: rank, comm_size 

  real*8    :: Fx, Fy, Fz

  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength
 
 
  !***********************************************************************
  !*                  intialisation (copied from jorek2_main)            *
  !***********************************************************************
#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif

  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
  
  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size
  
  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_GET_PROCESSOR_NAME (name,resultlength,ierr)
  write(*,'(A,I5,2A)') '  #MPI id, ProcessorName ', rank, ': ', name
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialize mode and mode_type arrays
  call det_modes()
 
  ! --- Preset input parameters to reasonable defaults, then read the input file.
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Initialize the vacuum part.
  call vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()

  if ( my_id == 0 ) then
    call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)
    if ( ierr /= 0 ) stop
  endif

  ! This is necessary for the parallel vacuum version during the code restart 
  call broadcast_phys(my_id)  
  if(freeboundary) call broadcast_vacuum(my_id, resistive_wall)

  call populate_element_rtree(node_list, element_list)
 
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Determine boundary information from the grid
  if ( my_id == 0 ) call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  
  ! --- Fill the vacuum response matrices for freeboundary computations
  if ( freeboundary ) then
    call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,    &
      resistive_wall)
    call update_response(my_id,tstep, freeboundary_equil, resistive_wall)
    call import_external_fields('coil_field.dat', my_id)
    call set_coil_curr_time_trace()
    call read_Z_axis_profile() 
    if ( (.not. restart) .or. (.not. wall_curr_initialized) ) call init_wall_currents(my_id, resistive_wall)
  end if
  
  call broadcast_elements(my_id, element_list)                ! elements

  call broadcast_nodes(my_id, node_list)                      ! nodes

  call populate_element_rtree(node_list, element_list)

  call broadcast_phys(my_id)                                  ! physics parameters

  if ( freeboundary ) call broadcast_vacuum(my_id, resistive_wall)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !***********************************************************************
  !*              end intialisation                                      *
  !***********************************************************************

  ! --- FORCES ---
  call total_wall_forces(my_id, node_list, element_list, Fx, Fy, Fz)

  if (my_id==0) then
    write(*,*) ' Fx = ', Fx
    write(*,*) ' Fy = ', Fy
    write(*,*) ' Fz = ', Fz
  endif

  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program jorek2_wall_forces 


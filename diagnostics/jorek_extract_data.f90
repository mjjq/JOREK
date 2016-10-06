!> Extract data from a restart file for non regression testing.
!!
!! Reads jorek_restart.rst and the namelist extract_data.nml. Writes an ascii file
!! extracted_data.dat which is compared to reference data for non-regression testing.
!!
!! See folder non-regression-testing/ in the repository.
!!
!! Supported values for extract_data:
!! * energies         -- Magnetic and thermal energies
!! * special_points   -- Axis and X-point positions and Psi values
!! * nodes            -- Geometry information from the node_list: node(i)%x(:,:)
program JOREK_EXTRACT_DATA
  
  use parameters
  use nodes_elements
  use phys_module
  use equil_info
  use mod_boundary
  use mod_log_params
  
  implicit none
  
  integer, parameter :: MAX_EXTRACT_DATA = 20
  
  integer                  :: ierr, i, j
  character(len=512)       :: extract_data(MAX_EXTRACT_DATA)
  namelist / extract / extract_data
  type (t_equil_state)     :: equil_state
  
  ! --- Initialization
  call det_modes()
  call initialise_parameters(0, "__NO_FILENAME__")
  call log_parameters(0)
  call initialise_basis
  call import_binary_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  
  ! --- If extract_data.nml file exists, read it.
  extract_data(:) = ''
  extract_data(1) = 'energies'
  open(42, file='./extract_data.nml', action='READ', status='OLD', iostat=ierr)
  if ( ierr == 0 ) then
    write(*,*) 'Reading parameters from extract_data.nml'
    read(42, extract)
    close(42)
  else
    write(*,*) 'extract_data.nml does not exist. Using default values.'
  end if
  
  ! --- Determine some information regarding the plasma state.
  call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, equil_state)
  call print_equil_state(equil_state, .true.)
  
  ! --- Write the requested data to the ascii file extracted_data.dat.
  open(42, file='./extracted_data.dat', action='WRITE', status='REPLACE')
  
  write(*,*) 'Extracting the following data:'
  
  do i = 1, MAX_EXTRACT_DATA
    
    if ( extract_data(i) == '' ) exit
    write(*,*) '  ', trim(adjustl(extract_data(i)))
    
    if ( extract_data(i) == 'energies' ) then
      
      write(42,'(es25.17)') energies(:,:,index_start)
      
    else if ( extract_data(i) == 'special_points' ) then
      
      write(42,'(es25.17)') equil_state%R_axis, equil_state%Z_axis, equil_state%Psi_axis
      if ( equil_state%xpoint ) then
        if ( ( equil_state%xcase == LOWER_XPOINT ) .or. ( equil_state%xcase == DOUBLE_NULL ) ) then
          write(42,'(es25.17)') equil_state%R_xpoint(1), equil_state%Z_xpoint(1),                  &
            equil_state%Psi_xpoint(1)
        else if ( ( equil_state%xcase == UPPER_XPOINT )                                            &
          .or. ( equil_state%xcase == DOUBLE_NULL ) ) then
          write(42,'(es25.17)') equil_state%R_xpoint(2), equil_state%Z_xpoint(2),                  &
            equil_state%Psi_xpoint(2)
        end if
      end if
      
    else if ( extract_data(i) == 'nodes' ) then
      
      do j = 1, node_list%n_nodes
        write(42,'(es25.17)') node_list%node(j)%x(:,:)
      end do
      
    else
      write(*,*) 'ERROR: Unknown value for extract_data.'
      stop
    end if
    
  end do
  
  close(42)
  
end program JOREK_EXTRACT_DATA

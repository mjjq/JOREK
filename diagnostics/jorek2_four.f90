PROGRAM JOREK2_FOUR
! For details see ../README.jorek2_four

  use phys_module
  use nodes_elements
  use fourier

  implicit none
  
  include 'mpif.h'

  type (type_surface_list) :: surface_list
  integer                  :: my_id
  integer                  :: i, j, k, l, ierr, ivar
  integer                  :: n_cpu, m_pol_range(2)
  integer                  :: required, provided, StatInfo
  integer*4                :: rank, comm_size 
  TYPE(t_theta_mapping)    :: mapping       ! mapping between theta_mag and theta_geo
  complex, allocatable     :: vfour(:,:,:,:)! Fourier transformed quantity.
  integer                  :: localvars(2), err, vars_per_cpu
  character(len=6)         :: sn
  ! Field line tracing parameters
  integer                  :: nstpts, nmaxsteps, nsmallsteps
  real                     :: deltaphi
  namelist / four_params / nstpts, nmaxsteps, deltaphi, nsmallsteps, m_pol_range

  required=MPI_THREAD_MULTIPLE
  call MPI_Init_thread(required,provided,StatInfo)         ! initialise threaded MPI (openMPI)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)           ! the id of each cpu
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)      ! the number of cpus
  my_id = rank
  n_cpu = comm_size
  
  !call system('date "+ %F %T.%N"')
  !WRITE(*,*) '@@> JOREK2_FOUR'
  
  call initialise_parameters(my_id)                        ! default values and namelist input
  call initialise_basis                                    ! define the basis functions at the Gaussian points

  call import_restart(node_list,element_list)              ! read restart file

  call broadcast_elements(my_id,element_list)              ! sending all elements
  call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)! sending boundary elements
  call broadcast_nodes(my_id,node_list)                    ! sending all nodes
  call broadcast_phys(my_id)                               ! sending the physics parameters
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Determine which variables shall be transformed by which MPI thread.
  vars_per_cpu = ( n_var - 1 ) / n_cpu + 1
  localvars = (/ my_id*vars_per_cpu+1, MIN( (my_id+1)*vars_per_cpu, n_var ) /)
    
  ! --- Preset field line tracing parameters.
  nstpts      = 30
  nmaxsteps   = 2500
  deltaphi    = 1.
  nsmallsteps = 3
  m_pol_range = (/1, 7/)
  
  ! --- If four_params.nml file exists, read field line tracing parameters from that file.
  OPEN(42, FILE='./four_params.nml', ACTION='READ', STATUS='OLD', IOSTAT=err)
  IF ( err == 0 ) THEN
    write(*,*) 'four_params.nml exists'
    READ(42, four_params) ! read namelist with Fourier parameters
    CLOSE(42)
  END IF
  
  ! --- Log field line tracing parameters.
  111 format(1x,a,2i12)
  112 format(1x,a,2es12.4)
  write(*,111) 'nstpts      =', nstpts
  write(*,111) 'nmaxsteps   =', nmaxsteps
  write(*,112) 'deltaphi    =', deltaphi
  write(*,111) 'nsmallsteps =', nsmallsteps
  write(*,111) 'm_pol_range =', m_pol_range

  ! --- Determine magnetic coordinates by field line tracing.
  CALL determine_theta_mag(nstpts, nmaxsteps, deltaphi, nsmallsteps, mapping, m_pol_range)
  
  ! --- Transform the quantities
  CALL transform_qttys(mapping, vfour, m_pol_range)
  
  ! --- Output Fourier modes of the physical quantities.
  do ivar = localvars(1), localvars(2)
    OPEN(43, FILE='mode_numbers', STATUS='REPLACE', ACTION='WRITE')
    write(42,'("# psi_normalized    ABS(",A," modes)")') TRIM(variable_names(ivar))
    write(42,'("#")')
    
    l = 0
    do j = 1, (n_tor+1)/2 ! tor
    
      write(sn,'(i3.3)') (j-1)*n_period
      OPEN(42, FILE=TRIM(variable_names(ivar))//'_modes_n'//trim(sn), STATUS='REPLACE', ACTION='WRITE')
    
      do i = m_pol_range(1), m_pol_range(2) ! pol
        
        write(42,'("# ",I3,":   m=",I3,", n=",I3)') l, i-1, (j-1)*n_period
        write(43,*) TRIM(int2str(i-1))//'/'//TRIM(int2str((j-1)*n_period))
        l = l + 1
        do k = 1, mapping.nstpts
          write(42,*) mapping.psin(k), ABS(vfour(i,j,k,ivar))
        end do
        write(42,*)
        write(42,*)
        
      end do
      
      CLOSE(42)
      
    end do
    
    CLOSE(43)
  end do

  !WRITE(*,*) '@@< JOREK2_FOUR'
  !call system('date "+ %F %T.%N"')

  call MPI_FINALIZE(IERR)                                ! clean up MPI
  
  
  
  
  
  
  CONTAINS
  
  
  
  
  
  
  CHARACTER(LEN=20) FUNCTION int2str(i)
    INTEGER, INTENT(IN) :: i
    
    write(int2str,*) i
    
    int2str = TRIM(ADJUSTL(int2str))
    
  END FUNCTION int2str
  
  
  
  
  
  
END PROGRAM JOREK2_FOUR

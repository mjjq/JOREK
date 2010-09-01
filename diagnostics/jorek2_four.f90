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
  integer                  :: n_cpu
  integer                  :: required, provided, StatInfo
  integer*4                :: rank, comm_size 
  TYPE(t_theta_mapping)    :: mapping       ! mapping between theta_mag and theta_geo
  real,    allocatable     :: vve(:,:,:)    ! Variable values at positions of mapping.rre and .zze
  complex, allocatable     :: vfour(:,:,:)  ! Fourier transformed quantity.
  integer                  :: localvars(2), err, vars_per_cpu
  ! Field line tracing parameters
  integer                  :: nstpts, nmaxsteps, nsmallsteps
  real                     :: deltaphi
  namelist / four_params / nstpts, nmaxsteps, deltaphi, nsmallsteps

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
  tstep = tstep_in

  call broadcast_elements(my_id,element_list)              ! sending all elements
  call broadcast_boundary(my_id,boundary_list)             ! sending boundary elements
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
  
  ! --- If four_params.nml file exists, read field line tracing parameters from that file.
  OPEN(42, FILE='./four_params.nml', ACTION='READ', STATUS='OLD', IOSTAT=err)
  IF ( err == 0 ) THEN
    write(*,*) 'four_params.nml exists'
    READ(42, four_params) ! read namelist with Fourier parameters
    CLOSE(42)
  END IF
  
  ! --- Log field line tracing parameters.
  write(*,*) 'nstpts      =', nstpts
  write(*,*) 'nmaxsteps   =', nmaxsteps
  write(*,*) 'deltaphi    =', deltaphi
  write(*,*) 'nsmallsteps =', nsmallsteps

  ! --- Determine magnetic coordinates by field line tracing.
  CALL determine_theta_mag(nstpts, nmaxsteps, deltaphi, nsmallsteps, mapping)

  ! --- Output Fourier modes of the physical quantities.
  do ivar = localvars(1), localvars(2)
    OPEN(42, FILE=TRIM(variable_names(ivar))//'_modes', STATUS='REPLACE', ACTION='WRITE')
    CALL transform_qtty(mapping, ivar, vve, vfour)
    
    l = 0
    do i = 1, 7 ! pol
      do j = 1, (n_tor+1)/2 ! tor
        
        write(42,'("# ",I3,":   m=",I3,", n=",I3," REAL")') l, i-1, j-1
        l = l + 1
        do k = 1, mapping.nstpts
          write(42,*) mapping.psin(k), REAL(vfour(i,j,k))
        end do
        write(42,*)
        write(42,*)
        
        write(42,'("# ",I3,":   m=",I3,", n=",I3," IMAG")') l, i-1, j-1
        l = l + 1
        do k = 1, mapping.nstpts
          write(42,*) mapping.psin(k), AIMAG(vfour(i,j,k))
        end do
        write(42,*)
        write(42,*)
        
      end do
    end do
    
    CLOSE(42)
  end do

  !WRITE(*,*) '@@< JOREK2_FOUR'
  !call system('date "+ %F %T.%N"')

  call MPI_FINALIZE(IERR)                                ! clean up MPI
END PROGRAM JOREK2_FOUR

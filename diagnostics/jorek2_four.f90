!> Allows to Fourier-analyse the physical variables of a JOREK restart file in magnetic coordinates.
program JOREK2_FOUR

  use parameters,     only: n_tor, n_var, n_period, variable_names
  use nodes_elements, only: element_list, node_list
  use fourier,        only: t_theta_mapping, determine_theta_mag, transform_qttys

  implicit none
  
  integer                  :: i, j, k, l, ierr, ivar, n_cpu, m_pol_range(2), err, vars_per_cpu
  type(t_theta_mapping)    :: mapping        ! Mapping between theta_mag and theta_geo
  complex, allocatable     :: vfour(:,:,:,:) ! Transformed quantities (m,n,irad,ivar)
  character(len=128)       :: filename
  
  ! ---Field line tracing parameters
  integer                  :: nstpts, nmaxsteps, nsmallsteps
  real                     :: deltaphi
  logical                  :: debug
  namelist / four_params / nstpts, nmaxsteps, deltaphi, nsmallsteps, m_pol_range, debug
  
  ! --- Initialization
  write(*,*)
  write(*,*) '>>> Initialization <<<'
  call initialise_parameters(0)                 ! default values and namelist input
  call initialise_basis                         ! define the basis functions at the Gaussian points
  call import_restart(node_list,element_list)   ! read restart file
  !   --- Preset field line tracing parameters.
  nstpts      = 30
  nmaxsteps   = 2500
  deltaphi    = 1.
  nsmallsteps = 3
  m_pol_range = (/0, 7/)
  debug       = .false.
  !   --- If four_params.nml file exists, read field line tracing parameters from that file.
  open(42, file='./four_params.nml', action='READ', status='OLD', iostat=err)
  if ( err == 0 ) then
    write(*,*) 'Reading parameters from four_params.nml'
    read(42, four_params) ! read namelist with Fourier parameters
    close(42)
  else
    write(*,*) 'WARNING: Could not find file four_params.nml -- using default parameters.'
  end if
  
  ! --- Log field line tracing parameters.
  write(*,*)
  write(*,*) '>>> Fourier transformation parameters <<<'
  111 format(1x,a,2i6)
  112 format(1x,a,2es12.4)
  113 format(1x,a,2l6)
  write(*,111) 'nstpts      =', nstpts
  write(*,111) 'nmaxsteps   =', nmaxsteps
  write(*,112) 'deltaphi    =', deltaphi
  write(*,111) 'nsmallsteps =', nsmallsteps
  write(*,111) 'm_pol_range =', m_pol_range
  write(*,113) 'debug       =', debug

  ! --- Determine magnetic coordinates by field line tracing.
  write(*,*)
  write(*,*) '>>> Determining the poloidal straight field line angle theta_star <<<'
  call determine_theta_mag(nstpts, nmaxsteps, deltaphi, nsmallsteps, mapping, m_pol_range, debug)
  
  ! --- Transform the quantities
  write(*,*)
  write(*,*) '>>> Transforming the variables <<<'
  call transform_qttys(mapping, vfour, m_pol_range)
  
  ! --- Output Fourier modes of the physical quantities.
  write(*,*)
  write(*,*) '>>> Writing the data to ascii files <<<'
  do ivar = 1, n_var
    do j = 1, (n_tor+1)/2 ! tor
    
      write(filename,'(a,a,i3.3)') trim(variable_names(ivar)), '_modes_n', (j-1)*n_period
      open(42, file=trim(filename), status='REPLACE', action='WRITE')
      l = 0
      do i = m_pol_range(1), m_pol_range(2) ! pol
        
        write(42,'("# ",I3,":   m=",I3,", n=",I3)') l, i-1, (j-1)*n_period
        l = l + 1
        do k = 1, mapping.nstpts
          write(42,*) mapping.psin(k), ABS(vfour(i,j,k,ivar))
        end do
        write(42,*)
        write(42,*)
        
      end do
      close(42)
      
    end do
  end do
  
  write(*,*) 'done.'
  
end program JOREK2_FOUR

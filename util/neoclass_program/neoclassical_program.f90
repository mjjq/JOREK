program neo

  use countlines_mod
  use prec_const
  use neo_subroutines

  implicit none

  integer 		:: nws, i
! INPUT READ FROM JOREK OUTPUT (CALCULATED AFTER EQUILIBRIUM)
  real(RKIND), dimension(:), allocatable        :: psi_pol
  real(RKIND), dimension(:), allocatable        :: T_read
  real(RKIND), dimension(:), allocatable        :: Ne_read
  real(RKIND), dimension(:), allocatable        :: q_read
! ---------------------------------------------------------                         
! USEFUL INPUT PARAMETERS
  real(RKIND)                                  :: central_density ! unit: m-3
  real(RKIND)                                  :: Bt
  real(RKIND)                                  :: R0
  real(RKIND)                                  :: amin
  character(512)                               :: qprofile_file
  character(512)                               :: T_Ne_profile_file
!----------------------------------------------------------
  real(RKIND)      	                       :: Qc
  real(RKIND)      	                       :: eps
  real(RKIND), dimension(:), allocatable       :: r
  real(RKIND), dimension(:), allocatable       :: psi_tor
  real(RKIND), dimension(:), allocatable       :: rho
  real(RKIND), dimension(:), allocatable       :: Ti
  real(RKIND), dimension(:), allocatable       :: Te
  real(RKIND), dimension(:), allocatable       :: Ne
  real(RKIND), dimension(:), allocatable       :: Ni
  real(RKIND), dimension(:), allocatable       :: q
! real(RKIND), dimension(:), allocatable       :: ke
  real(RKIND), dimension(:), allocatable       :: ki
! real(RKIND), dimension(:), allocatable       :: muneo_e
  real(RKIND), dimension(:), allocatable       :: muneo_i
  real(RKIND), dimension(:), allocatable       :: muneo_i_jorek
  real(RKIND), dimension(:,:), allocatable     :: muneoe
  real(RKIND), dimension(:,:), allocatable     :: muneoi
  integer                                      :: err



! READ INPUT PARAMETERS: q, T, Ne, psi_pol

  ! --- Namelist with input parameters.
  namelist /in1/ central_density, Bt, R0, amin,     &
                 qprofile_file, T_Ne_profile_file

  read(5,in1)

  OPEN(UNIT=72, FILE=trim(qprofile_file), FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
  if ( err /= 0 ) then
     write(*,*) 'ERROR: Cannot open file '//TRIM(qprofile_file)//'.'
     stop
  endif
  OPEN(UNIT=73, FILE=trim(T_Ne_profile_file), FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
  if ( err /= 0 ) then
     write(*,*) 'ERROR: Cannot open file '//TRIM(T_Ne_profile_file)//'.'
     stop
  endif

  nws = countlines(72)
  print*, 'nws = ', nws
  if (nws .lt. 2) then
     write(*,*) '  ERROR: Could not read the numerical profiles ''qprofile'' '
     stop
  end if
  i = countlines(73)
  if (i .ne. nws) then
     write(*,*) '  ERROR: the profiles do not have the same size'
     stop
  endif

  allocate(psi_pol(nws), T_read(nws), Ne_read(nws), q_read(nws))
  do i=1,nws
     read(72,*) psi_pol(i), q_read(i)
     read(73,*) T_read(i), Ne_read(i)
  end do

  close(72)
  close(73)

  !read(74,*) central_density, Bt, R0, amin
!!$  central_density = 1.0d20
!!$  Bt = 5.3d0
!!$  R0 = 6.19476d0
!!$  amin = 2.d0

!  allocate(r(nws), Te(nws), Ti(nws), Ne(nws), Ni(nws), rho(nws), muneoe(nws,3), muneoi(nws,3))
  allocate(r(nws), Te(nws), Ne(nws), psi_tor(nws), rho(nws))
  allocate(muneoe(nws,3), muneoi(nws,3))
  allocate(muneo_i(nws), muneo_i_jorek(nws), ki(nws))
  r=sqrt(psi_pol)

! Calculation of the toroidal flux psi_tor (we need rho=sqrt(psi_tor) to calculate the factor of trapped particles in calculs_coeffs.f90)
  psi_tor(1)=0
  DO i=1,nws-1
     Qc=0.5*(q_read(i)+q_read(i+1))
     psi_tor(i+1)=psi_tor(i)+Qc*((r(i+1))**2-(r(i))**2)
  END DO
  
  DO i=1,nws
     rho(i)=sqrt(psi_tor(i)/psi_tor(nws))
  END DO
  
! de-normalization from JOREK
  Ne = Ne_read* central_density * 1.d20 !!! central_density in 1.d20 m-3
!  Ni = Ne
  Te = T_read / (2 * charge * mu0 * central_density * 1.d20)
!  Ti = Te

  eps = R0/amin

! calculation of neoclassical coefficients depending on r=sqrt(psi_pol)
!  call neocoeffs(nws, eps, Te, Ti, Ne, Ni, q_read, rho, muneoe, muneoi)
  call neocoeffs(nws, eps, R0, Te, Te, Ne, Ne, q_read, rho, muneoe, muneoi)

  do i=1,nws
     muneo_i(i)=muneoi(i,1)
     ki(i)=muneoi(i,2)/muneoi(i,1)
  end do

  !normalization of muneo_i
  muneo_i_jorek = muneo_i * sqrt(mu0*mion*central_density*1.d20)
  muneo_i_jorek(1) = 0.d0
  ki(1)=ki(2)

! create a new file with muneo_i_jorek and ki radial profiles
  OPEN(75,FILE='neoclass_coef.dat', status='replace', action='write')
  write (75,*)    0.d0,        0.d0,         0.d0
  do i=1,nws
     write (75,*) psi_pol(i), muneo_i_jorek(i), ki(i)
  end do
  write(*,*) 'file ''neoclass_coef.dat'' created containing psi_pol, mu_neo_i and k_i'

  deallocate(psi_pol, q_read, T_read, Ne_read, r, rho, psi_tor, Te, Ne, muneoi, muneoe, muneo_i, muneo_i_jorek, ki)    

STOP
end program neo

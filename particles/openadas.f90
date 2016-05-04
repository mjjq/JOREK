!------------------------------------------------------------------------------------------------------
! module takes the OPEN-ADAS data to calculate the steady state (or time evolution) charge distribution
! and average charge state as a function of temperature
!------------------------------------------------------------------------------------------------------

module openadas

type type_ADF11 !< Custom data structure containing relevant fields from ADF11 format files (unresolved case!)
  integer             :: n_Z !< Atomic number
  integer             :: izmin, izmax !< minimum and maximum value of z for which data is available
  real*8, allocatable :: density(:) !< log10 density (m^-3)
  real*8, allocatable :: temperature(:) !< log10 temperature (K)
  real*8, allocatable :: GRC(:,:,:) !< log10 of coefficient (parameters: z, d, T). Units:
  ! ACD, SCD: cm3s-1 (for *CD ?)
  ! PLT, PRB: Wcm3 (for P* ?)
end type type_ADF11

type type_ADF11_all ! unresolved
  integer          :: n_Z !< Atomic number
  type(type_ADF11) :: ACD !< Effective recombination coefficients
  type(type_ADF11) :: SCD !< Effective ionisation coefficients
  type(type_ADF11) :: CCD !< Charge exchange effective recombination coefficients
  type(type_ADF11) :: PLT !< Line power driven by excitation of dominant ions
  type(type_ADF11) :: PRB !< Continuum and line power driven by recombination and bremsstrahlung of dominant ions
  type(type_ADF11) :: PRC !< Line power due to charge transfer from thermal neutral hydrogen to dominant ions
end type type_ADF11_all
!! Recombination data is given as recombining FROM (Z=1 to Z=74)
!! Ionisation data is given as ionising TO (Z=1 up to Z=74)
contains

!> Read ADF11 data files and import them into a type_ADF11
!! Tries to read ACD, SCD, CCD, PLT, PRB, PRC coefficients
!! if the files exist. Files of format acd$suffix.dat are read.
!! Suffix is usually of the form 50_w, 96_li
function read_adf11(suffix) result(ad)
use constants
implicit none

character*6, intent(in) :: suffix !< Usually year_atom (ex: 50_w, 96_li), give this in input file
type(type_ADF11_all), target :: ad !< OpenAdas data type

type(type_ADF11), pointer :: a
integer :: i_ADF11
character*3, dimension(1:6), parameter :: ADF11_filenames = (/"acd", "scd", "ccd", "plt", "prb", "prc"/)
character*13 :: filename

integer :: i, ierr, n_d, n_T
logical :: file_exists

write(*,'(A)') '*********************************'
write(*,'(A)') '* Importing OpenAdas data       *'
write(*,'(A,A,A)') '* open files ending in: ', suffix, '  *' 
write(*,'(A)') '*********************************'

do i_ADF11 = 1,size(ADF11_filenames,1)
  write(filename,"(A,A,A)") ADF11_filenames(i_ADF11), trim(suffix), '.dat'
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) cycle ! Skip this type of data

  write(*,"(A,A)",advance="no") "Reading data from ", filename
  open(10,file=filename,status="old",iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) "failed with code", ierr
    cycle
  endif

  ! Point a to right variable to read in
  select case (i_ADF11)
    case (1); a => ad%ACD
    case (2); a => ad%SCD
    case (3); a => ad%CCD; write(*,*) "Warning: CCD not implemented correctly yet"
    case (4); a => ad%PLT
    case (5); a => ad%PRB
    case (6); a => ad%PRC; write(*,*) "Warning: PRC not implemented correctly yet" ! see coronal model
  end select

  read(10,*)  a%n_z, n_d, n_T, a%izmin, a%izmax
  ad%n_z = a%n_z
  allocate(a%density(n_d), a%temperature(n_T), a%GRC(a%izmin:a%izmax,n_d,n_T))

  read(10,*)
  read(10,*) a%density(:)
  ! Convert densities to log10 of m^-3 instead of cm^-3
  a%density = a%density + 6.d0
  read(10,*) a%temperature(:)
  ! Convert temperatures to log10 of K instead of eV
  ! From E eV = kB T
  a%temperature = a%temperature - log10(K_BOLTZ) + log10(EL_CHG) ! EL_CHG * 1 Volt actually

  a%GRC = -30.d0
  do i = a%izmin, a%izmax
    read(10,*)
    read(10,*) a%GRC(i,:,:)
  enddo
  close(10)

  write(*,"(A)") "succeeded"
enddo

! Test if ACD and SCD were loaded at least
if (.not. (allocated(ad%ACD%density) .and. allocated(ad%SCD%density))) then
  write(*,*) "ACD and SCD not found, exiting"
  call exit(10)
else
  write(*,*) 'Done reading adas data for atomic number', ad%n_Z
endif
end function read_adf11



!> (Zeroth order) interpolation of log10 values of GRC in density and temperature
!! (Just return closest value for now, replace with splines in the future)
pure function GRC(a, z, density, temperature)
implicit none

type (type_ADF11), intent(in) :: a
real*8, intent(in)            :: density     !< log10 density in m^-3
real*8, intent(in)            :: temperature !< log10 temperature in K
integer, intent(in)           :: z !< index in a%GRC(z,:,:) (is ionisation level or ionisation level - 1, 1:n_z)
real*8 :: GRC

integer :: index_d, index_T

if (allocated(a%GRC)) then
  index_d = minloc(abs(density-a%density),1)
  index_T = minloc(abs(temperature-a%temperature),1)

  GRC = 10.d0**a%GRC(z,index_d,index_T)
else
  GRC = 0.d0
endif
end function GRC


end module openadas

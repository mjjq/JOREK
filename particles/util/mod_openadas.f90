!> module takes the OPEN-ADAS data to calculate the steady state (or time evolution) charge distribution
!> and average charge state as a function of temperature
module mod_openadas

!> Custom data structure containing relevant fields from ADF11 format files (unresolved case!)
type type_ADF11
  integer             :: n_Z !< Atomic number
  integer             :: izmin, izmax !< minimum and maximum value of z for which data is available
  real*8, allocatable :: density(:) !< log10 density (m^-3)
  real*8, allocatable :: temperature(:) !< log10 temperature (K)
  real*8, allocatable :: GRC(:,:,:) !< log10 of coefficient (parameters: d, T, z). Units:
  !< ACD, SCD: cm3s-1 (for *CD ?)
  !< PLT, PRB: Wcm3 (for P* ?)
end type type_ADF11

!> Compound datatype containing many type_ADF11
type type_ADF11_all
  integer          :: n_Z !< Atomic number
  type(type_ADF11) :: ACD !< Effective recombination coefficients
  type(type_ADF11) :: SCD !< Effective ionisation coefficients
  type(type_ADF11) :: CCD !< Charge exchange effective recombination coefficients
  type(type_ADF11) :: PLT !< Line power driven by excitation of dominant ions
  type(type_ADF11) :: PRB !< Continuum and line power driven by recombination and bremsstrahlung of dominant ions
  type(type_ADF11) :: PRC !< Line power due to charge transfer from thermal neutral hydrogen to dominant ions
end type type_ADF11_all
!< Recombination data is given as recombining FROM (Z=1 to Z=74)
!< Ionisation data is given as ionising TO (Z=1 up to Z=74)
contains

!> Read ADF11 data files and import them into a type_ADF11
!> Tries to read ACD, SCD, CCD, PLT, PRB, PRC coefficients
!> if the files exist. Files of format acd$suffix.dat are read.
!> Suffix is usually of the form 50_w, 96_li
function read_adf11(suffix) result(ad)
use constants
use_mpi
implicit none

character*6, intent(in) :: suffix !< Usually year_atom (ex: 50_w, 96_li), give this in input file
type(type_ADF11_all), target :: ad !< OpenAdas data type

type(type_ADF11), pointer :: a
integer :: i_ADF11
character*3, dimension(1:6), parameter :: ADF11_filenames = (/"acd", "scd", "ccd", "plt", "prb", "prc"/)
character*13 :: filename

integer :: i, ierr, n_d, n_T, k, my_id
logical :: file_exists

call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

if (my_id .eq. 0) then
  write(*,'(A)') '*********************************'
  write(*,'(A)') '* Importing OpenAdas data       *'
  write(*,'(A,A,A)') '* open files ending in: ', suffix, '  *'
  write(*,'(A)') '*********************************'
endif

do i_ADF11 = 1,size(ADF11_filenames,1)
  write(filename,"(A,A,A)") ADF11_filenames(i_ADF11), trim(suffix), '.dat'
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) cycle ! Skip this type of data

  if (my_id .eq. 0) write(*,"(A,A)",advance="no") "Reading data from ", filename
  open(10,file=filename,status="old",iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) my_id, " failed with code ", ierr
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
  allocate(a%density(n_d), a%temperature(n_T), a%GRC(n_d,n_T,a%n_z))

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
    read(10,*) a%GRC(:,:,i)
  enddo
  close(10)

  if (my_id .eq. 0) write(*,"(A)") "succeeded"
enddo

! Test if ACD and SCD were loaded at least
if (.not. (allocated(ad%ACD%density) .and. allocated(ad%SCD%density))) then
  write(*,*) my_id, "ACD and SCD not found, exiting"
  call exit(10)
else
  if (my_id .eq. 0) write(*,*) 'Done reading adas data for atomic number', ad%n_Z
endif
end function read_adf11



!> interpolation of log10 values of GRC in density and temperature
pure function GRC(a, z, density, temperature)
implicit none

type (type_ADF11), intent(in) :: a           !< ADF11 datatype
real*8, intent(in)            :: density     !< log10 density in m^-3
real*8, intent(in)            :: temperature !< log10 temperature in K
integer, intent(in)           :: z !< index in a%GRC(:,:,z) (is ionisation level or ionisation level - 1, 1:n_z)
real*8 :: GRC !< Generalized Radiational Coefficient at this density and temperature

! If GRC exists and we are looking for a Z that is nonzero
if (allocated(a%GRC) .and. z .le. ubound(a%GRC,3) .and. z .ge. lbound(a%GRC,3)) then
  GRC = 10.d0**L2Dinterp(a%density,a%temperature,a%GRC(:,:,z),density,temperature)
else
  GRC = 0.d0
endif
end function GRC

!> Linear 2D interpolation on a rectangular grid
!> x2y1       xy1    x1y1
!>  *----------*------*
!>             |
!>             * xy
!>             |
!>             |
!>  *----------*------*
!> x2y2       xy2    x1y2
!>
!> Calculates the interpolation using two intermediate values
!> fxy1 and fxy2.
!> Equations used are:
!> \[fx1  = \frac{f_{11}-f_{21}}{x_1-x_2} (x-x_1) + f_{11}\]
!> \[fx2  = \frac{f_{12}-f_{22}}{x_1-x_2} (x-x_1) + f_{12}\]
!> \[fout = \frac{f_{x1}-f_{x2}}{y_1-y_2} (y-y_1) + f_{x1}\]
!> x1,2 and y1,2 are chosen in order of closeness
!> This algorithm can also be used for extrapolation
pure function L2Dinterp(tx,ty,f,x,y) result(fout)
implicit none

real*8, intent(in), dimension(:)                 :: tx !< Grid points in x
real*8, intent(in), dimension(:)                 :: ty !< Grid points in y
real*8, intent(in), dimension(size(tx),size(ty)) :: f !< Function values at these points
real*8, intent(in)  :: x, y !< Points at which to interpolate
real*8              :: fout

integer :: ix1, iy1 !< Index of closest point
integer :: ix2, iy2 !< Index of other (usually next closest) point
real*8  :: fx1, fx2 ! Temporary variables
ix1 = minloc(abs(tx - x), dim=1)
if (x .ge. tx(ix1)) ix2 = ix1 + 1 ! find other index
if (x .lt. tx(ix1)) ix2 = ix1 - 1
if (ix2 .gt. size(tx)) ix2 = size(tx) - 1 ! if it does not exist, extrapolate
if (ix2 .lt. 1       ) ix2 = 2
iy1 = minloc(abs(ty - y), dim=1)
if (y .ge. ty(iy1)) iy2 = iy1 + 1
if (y .lt. ty(iy1)) iy2 = iy1 - 1
if (iy2 .gt. size(ty)) iy2 = size(ty) - 1
if (iy2 .lt. 1       ) iy2 = 2

fx1  = (f(ix1,iy1) - f(ix2,iy1))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy1)
fx2  = (f(ix1,iy2) - f(ix2,iy2))/(tx(ix1) - tx(ix2)) * (x - tx(ix1)) + f(ix1,iy2)
fout = (fx1 - fx2) / (ty(iy1) - ty(iy2)) * (y - ty(iy1)) + fx1
end function L2Dinterp
end module mod_openadas

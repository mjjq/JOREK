!------------------------------------------------------------------------------------------------------
! module takes the OPEN-ADAS data to calculate the steady state (or time evolution) charge distribution
! and average charge state as a function of temperature
!------------------------------------------------------------------------------------------------------

module openadas

type type_ADF11 !< Custom data structure containing relevant fields from ADF11 format files (unresolved case!)
  integer             :: n_Z !< Atomic number
  integer             :: izmin, izmax !< minimum and maximum value of z for which data is available
  real*8, allocatable :: density(:) !< log10 density (cm^-3)
  real*8, allocatable :: temperature(:) !< log10 temperature (eV)
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

type type_coronal
  integer :: n_Z !< Atomic number
  real*8, allocatable :: density(:) !< log10 density (m^-3)
  real*8, allocatable :: temperature(:) !< log10 temperature (K)
  real*8, allocatable :: Z(:,:) !< Charge state (electroncharges)
  real*8, allocatable :: Prad(:,:) !< log10 Radiated power per ion (W)
end type type_coronal

contains

!> Read ADF11 data files and import them into a type_ADF11
!! Tries to read ACD, SCD, CCD, PLT, PRB, PRC coefficients
!! if the files exist. Files of format acd$suffix.dat are read.
!! Suffix is usually of the form 50_w, 96_li
function read_adf11(suffix) result(ad)
implicit none

character*6, intent(in) :: suffix !< Usually year_atom (ex: 50_w, 96_li), give this in input file
type(type_ADF11_all) :: ad !< OpenAdas data type

type(type_ADF11) :: a
integer :: i_ADF11
character*3, dimension(1:6), parameter :: ADF11_filenames = (/"acd", "scd", "ccd", "plt", "prb", "prc"/)
character*13 :: filename

integer :: i, ierr, n_d, n_T
logical :: file_exists

write(*,'(A)') '*********************************'
write(*,'(A)') '* Importing OpenAdas data       *'
write(*,'(A)') '* open files ending in:', suffix, ' *' 
write(*,'(A)') '*********************************'

do i_ADF11 = 1,size(ADF11_filenames,1)
  write(filename,"(A,A,A)") ADF11_filenames(i_ADF11), trim(suffix), '.dat'
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) continue ! Skip this type of data

  write(*,"(A,A,A)",advance="no") "Reading data from ", filename
  open(10,file=filename,iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) "failed with code", ierr
    continue
  endif

  read(10,*)  a%n_z, n_d, n_T, a%izmin, a%izmax
  allocate(a%density(n_d), a%temperature(n_T), a%GRC(a%izmin:a%izmax,n_d,n_T))

  read(10,*)
  read(10,*) a%density(:)
  read(10,*) a%temperature(:)

  a%GRC = -30.d0
  do i = a%izmin, a%izmax
    read(10,*)
    read(10,*) a%GRC(i,:,:)
  enddo
  close(10)

  select case (i_ADF11)
    case (1); ad%ACD = a
    case (2); ad%SCD = a
    case (3); ad%CCD = a; write(*,*) "Warning: CCD not implemented correctly yet"
    case (4); ad%PLT = a
    case (5); ad%PRB = a
    case (6); ad%PRC = a; write(*,*) "Warning: PRC not implemented correctly yet" ! see coronal model
  end select

  write(*,*) "succeeded"
enddo

write(*,*) 'Done reading adas data for atomic number', ad%n_Z
end function read_adf11



!> (Zeroth order) interpolation of log10 values of GRC in density and temperature
!! (Just return closest value for now, replace with splines in the future)
function GRC(a, z, density, temperature)
implicit none

type (type_ADF11), intent(in) :: a
real*8, intent(in)            :: density, temperature
integer, intent(in)           :: z ! index in a%GRC(z,:,:) (is ionisation level or ionisation level - 1, 1:n_z)
real*8 :: GRC

integer :: index_d, index_T

index_d = minloc(abs(density-a%density),1)
index_T = minloc(abs(temperature-a%temperature),1)

if (allocated(a%GRC)) then
  GRC = 10.d0**a%GRC(z,index_d,index_T)
else
  GRC = 0.d0
endif
end function GRC



!> Calculate the coronal equilibrium values at specific values of density and temperature
!! 
function coronal_equilibrium(ad) result(cor)
implicit none

type (type_ADF11_all), intent(in) :: ad
type (type_coronal)               :: cor

real*8, allocatable :: coronal_Z(:), corona_matrix(:,:), zni(:), b(:), Z(:), fractions(:)
real*8, allocatable :: A_l(:),A_d(:),A_u(:),T_bar(:),Z_bar(:), R_bar(:)
real*8              :: density, temperature, log_density, log_temperature

real*8, allocatable :: ion_rate(:), rec_rate(:), rad_rate(:), Z_states(:)
real*8              :: zc(30)
integer             :: n_states, i, info, n_step, j, k, m
integer             :: i_d_ion, i_T_ion, i_d_rec, i_T_rec

integer :: n_d, n_T, iz

n_states = ad%n_Z + 1

write(*,'(A)')      '*********************************'
write(*,'(A,i3,A)') '* Coronal model : ',n_states,'           *'
write(*,'(A)')      '*********************************'

cor%n_Z = ad%n_Z
n_d = 4
n_T = 200

allocate(cor%density(n_d), cor%temperature(n_T), cor%Z(n_d,n_T), cor%Prad(n_d,n_T))

allocate(corona_matrix(3,n_states),zni(n_states),b(n_states), Z(n_states), fractions(n_states))

allocate(A_l(n_states-1),A_d(n_states),A_u(n_states-1))

cor%density = (/ 18., 19., 20., 21. /) ! log10 [m^-3]

do m=1, n_d

  density = 10**cor%density(m) ! m^-3

  do k=1, n_T

    cor%temperature(k) = log10( 1.d0 + exp(log(4.d4)*float(k-1)/(float(n_T-1))) - 1.d0 )  ! in [K]

    do iz=1,cor%n_Z
      ion_rate(iz) = GRC(ad%SCD, iz, cor%density(m), cor%temperature(k)) ! ionizing to level iz
      rec_rate(iz) = GRC(ad%ACD, iz, cor%density(m), cor%temperature(k)) + &
                     GRC(ad%CCD, iz, cor%density(m), cor%temperature(k)) ! recombining from level iz
      ! These should actually be multiplied by the electron and hydrogen ion densities
      ! In the trace impurity limit these are the same and can be skipped (multiplication by constant does not change fractions)
      rad_rate(iz) = GRC(ad%PRB, iz, cor%density(m), cor%temperature(k)) + &
                     GRC(ad%PLT, iz, cor%density(m), cor%temperature(k)) + &
                     GRC(ad%PRC, iz, cor%density(m), cor%temperature(k)) ! radiation emitted by atoms at level iz
      ! PRB and PLT should also be multiplied by n_e, and PRC with neutral density
    end do

    corona_matrix = 0.d0

    ! Create diagonal elements of coronal model matrix
    do i=2, cor%n_Z
      Z(i) = float(i-1)
      corona_matrix(1,i) = + ion_rate(i-1) ! increase from ionisation
      corona_matrix(2,i) = - ion_rate(i)   - rec_rate(i) ! Loss in this state
      corona_matrix(3,i) =                 + rec_rate(i+1) ! increase from recombination of higher level atoms
    enddo
    Z(1)        = 0.
    Z(cor%n_Z+1) = float(cor%n_Z)

    corona_matrix(1,1)     = 0.
    corona_matrix(2,1)     = - ion_rate(1)      ! ionisation losses from n=0 to 1
    corona_matrix(3,1)     = + rec_rate(1)      ! recombination from n=1 to 0

    corona_matrix(1,cor%n_Z+1) = + ion_rate(cor%n_Z)
    corona_matrix(2,cor%n_Z+1) = - rec_rate(cor%n_Z)
    corona_matrix(3,cor%n_Z+1) = 0.

    corona_matrix = corona_matrix * density

    zni = 1.e9

    do i=2,cor%n_Z
      b(i) =  zni(i-1)*corona_matrix(1,i) + zni(i)*corona_matrix(2,i) + zni(i+1)*corona_matrix(3,i)
    enddo
    b(1)        = zni(1)*corona_matrix(2,1)                 + zni(2)*corona_matrix(3,1)
    b(n_states) = zni(n_states-1)*corona_matrix(1,n_states) + zni(n_states)*corona_matrix(2,n_states)

    ! Lower, diagonal and upper components
    A_l  =    - corona_matrix(1,2:n_states)
    A_d  = 1. - corona_matrix(2,:)
    A_u  =    - corona_matrix(3,1:n_states-1)

    ! Right hand side
    b = b

    ! Solve AX=B for tridiagonal matrices. Result stored in b
    call dgtsv(n_states,1,A_l,A_d,A_u,b,n_states,info)

    zni = zni + b ! minimum levels of 1e9?

    if (info .ne. 0) write(*,*) 'info : ',info

    fractions    = zni / sum(zni)
    cor%Z(m,k)   = dot_product(fractions,Z)
    cor%Prad(m,k) = log10( dot_product(fractions,rad_rate))
  enddo
enddo
end function coronal_equilibrium



!> (Zeroth order) interpolation of coronal model charge
!! (Just return closest value for now, replace with splines in the future)
subroutine interpolate_coronal(cor, density, temperature, z, rad)
implicit none

type (type_coronal), intent(in) :: cor
real*8, intent(in)              :: density ! log10 density (m^-3)
real*8, intent(in)              :: temperature ! log10 temperature (K)
real*8, intent(out)             :: z ! most probable charge state
real*8, intent(out)             :: rad ! radiated power according to coronal equilibrium

integer :: index_d, index_T

index_d = minloc(abs(density-cor%density),1)
index_T = minloc(abs(temperature-cor%temperature),1)

z   = cor%Z(index_d,index_T)
rad = 10.d0**cor%Prad(index_d,index_T) ! Still needs to be multiplied by electron density! XXX
end subroutine interpolate_coronal

end module openadas

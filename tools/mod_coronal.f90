!> module takes the OPEN-ADAS data to calculate the coronal equil-
!> ibrium temperature and radiation, optionally in time.
module mod_coronal
use mod_openadas
implicit none
private
public coronal
public coronal_timestep ! tested in [[compare_mc_coronal]]

!> Coronal equilibrium datatype
type coronal
  integer :: n_Z !< Atomic number
  real*8, allocatable :: density(:) !< log10 density (m^-3)
  real*8, allocatable :: temperature(:) !< log10 temperature (K)
  real*8, allocatable :: Z(:,:,:) !< Charge state (e) density for specific densities, temperatures and charge states [i_n, i_T, i_q]
  real*8, allocatable :: Prad(:,:) !< log10 Radiated power per ion (W) for the above densities and temperatures [i_n, i_T]
contains
  procedure :: interp => interpolate_coronal
  procedure :: interp_gradients => interpolate_coronal_gradients
end type coronal
interface coronal
  module procedure coronal_equilibrium
end interface coronal

contains

!> Calculate the coronal matrix in tridiagonal form
!> This matrix defines the time derivatives of population sizes
!> As drho/dt = A rho where A is this matrix.
!> See http://jorek.eu/wiki/doku.php?id=adf11 for more info
pure function corona_matrix(ad, density, temperature)
type (ADF11_all), intent(in)         :: ad !< ADF11 datatype
real*8, intent(in)                   :: density !< log10 density in m^-3
real*8, intent(in)                   :: temperature !< log10 temperature in K
real*8, dimension(1:3,0:ad%n_Z)      :: corona_matrix !< Output the coronal matrix diagonals

real*8, dimension(0:ad%n_Z+1) :: ion_rate, rec_rate
integer :: iz, i

! Get the ionization and recombination
ion_rate = 0.d0
rec_rate = 0.d0
do iz=1,ad%n_Z
  ion_rate(iz) = ad%SCD%interp(iz, density, temperature) ! ionizing to level iz (0 is neutral)
  rec_rate(iz) = ad%ACD%interp(iz, density, temperature) + &
                 ad%CCD%interp(iz, density, temperature) ! recombining from level iz
  ! These should actually be multiplied by the electron and hydrogen ion densities (done below)
  ! In the trace impurity limit these are the same
end do

! Fill in tridiagonal elements (implicitly shifted down one row per column)
do i=0,ad%n_Z
  corona_matrix(1,i) = + ion_rate(i) ! increase from ionisation
  corona_matrix(2,i) = - ion_rate(i+1) - rec_rate(i) ! Loss in this state
  corona_matrix(3,i) =                 + rec_rate(i+1) ! increase from recombination of higher level atoms
enddo

corona_matrix = corona_matrix * 10.d0**density
end function corona_matrix


!> Radiated power in a specific coronal equilibrium configuration and temperature
pure function coronal_Prad(ad, density, temperature, fractions, neutral_density)
type (ADF11_all), intent(in)            :: ad !< ADF11 datatype
real*8, intent(in)                      :: density !< log10 density in m^-3
real*8, intent(in)                      :: temperature !< log10 electron temperature in K
real*8, intent(in), dimension(0:ad%n_Z) :: fractions !< Fractional charge states. Should sum to 1 but we do not check it!
real*8, intent(in), optional            :: neutral_density !< log10 neutral density in m^-3
real*8                                  :: coronal_Prad !< Output power in W / atom

real*8, dimension(ad%n_Z) :: rad
real*8, dimension(ad%n_Z) :: rad_RC
real*8  :: density_n
integer :: iz

do iz=1,ad%n_Z
  rad(iz)      = ad%PRB%interp(iz, density, temperature) + &
                 ad%PLT%interp(iz, density, temperature)
  rad_RC(iz)   = ad%PRC%interp(iz, density, temperature)
enddo ! radiation emitted by atoms at level iz
! PRB and PLT should also be multiplied by n_e, and PRC with neutral density

! If the neutral Hydrogen density is present use it
if (present(neutral_density)) then
  density_n = neutral_density
else ! otherwise set it to some extremely low value
  density_n = -99.d0 ! this is log10 of density
endif

coronal_Prad = dot_product(fractions(1:ad%n_Z), rad*10.d0**density + rad_RC*10.d0**density_n)
end function coronal_Prad


!> Calculate the coronal equilibrium values at specific values of density and temperature
function coronal_equilibrium(ad) result(cor)
use constants
type (ADF11_all), intent(in) :: ad !< ADF11 datatype
type (coronal)               :: cor !< Coronal equilibrium datatype

real*8, dimension(0:ad%n_Z) :: p, p_old, Z
integer :: n_d, n_T, iz, k, m, i, j
real*8, parameter :: tolerance = 1d-12 ! on max(abs(new - old))
logical :: converged
real*8, dimension(1:3,0:ad%n_Z) :: cmat
real*8 :: dt

cor%n_Z = ad%n_Z
n_d = 10
n_T = 200
do iz=0,ad%n_Z
  Z(iz) = real(iz,8)
enddo

allocate(cor%density(n_d), cor%temperature(n_T), cor%Z(n_d,n_T,0:cor%n_Z), cor%Prad(n_d,n_T))
!$omp parallel do default(none) shared(cor, ad, n_d, n_T, Z) private(p, p_old, m, k, converged, cmat, i, j, dt)
do m=1, n_d
  cor%density(m) = 18.d0 + float(m-1)/n_d * (21.-18.) ! log10 [m^-3], linear between 18 and 21

  ! Here we are using the distribution from the previous temperature as an
  ! initial guess of the next temperature. So only initialize at a different
  ! density scan instead of every temperature scan.
  p = 0.d0
  p(0) = 1.d0

  do k=1, n_T
    cor%temperature(k) = log10( 1.d0 + exp(log(4.d4)*float(k-1)/(float(n_T-1))) - 1.d0 ) + log10(EL_CHG) - log10(K_BOLTZ) ! in log10 [K]

    converged = .false.
    ! Iterate until converged or we run out of steps
    cmat = corona_matrix(ad, cor%density(m), cor%temperature(k))
    pow_loop: do i=0,5
      dt = 10.d0**(-i)
      step_loop: do j=1,10**i + 100 ! extra 100 to perform more large steps
        p_old = p
        call coronal_timestep(ad, p, dt, cmat) ! use a fixed timestep of 0.1s to solve the equilibrium state
        if (maxval(abs(p_old - p)) .lt. tolerance) then
          converged = .true.
          exit pow_loop
        end if
      end do step_loop
    end do pow_loop

    !$omp critical
    if (.not. converged) write(*,*) "WARNING: coronal equilibrium not converged, max delta", maxval(abs(p_old - p))
    !$omp end critical

    cor%Z(m,k,:)  = p/sum(p)
    cor%Prad(m,k) = coronal_Prad(ad, cor%density(m), cor%temperature(k), p/sum(p)) ! Do not set neutral density yet
  enddo
enddo
!$omp end parallel do
end function coronal_equilibrium


!> Perform one timestep with the coronal equilibrium matrix
!>
!> Equation to solve is: \(p' = A p\)
!> Discretize as \(p^{n+1} - p^n = \Delta t \left((1-\theta)Ap^n + \theta A p^{n+1}\right)\)
!> Leading to a matrix equation
!> \((1-\theta \Delta t A)p^{n+1}a = (1 + \Delta t (1-\theta) A) p^n \)
subroutine coronal_timestep(ad, p, tstep, cmat)
type (ADF11_all), intent(in)               :: ad !< ADF11 datatype
real*8, dimension(0:ad%n_Z), intent(inout) :: p !< Population in each level
real*8, intent(in)                         :: tstep !< Timestep size in s
real*8, dimension(1:3,0:ad%n_Z)            :: cmat !< tridiagonal corona matrix

real*8, parameter :: theta = 1.d0 ! 0.5 = Crank-Nicholson, 1.d0 = Backward Euler

real*8, dimension(0:ad%n_Z)     :: b
real*8 :: A_l(1:ad%n_Z), &
          A_d(0:ad%n_Z), &
          A_u(1:ad%n_Z)
integer :: info, i

! Lower, diagonal and upper components (I - theta dt A)
A_l =      - theta * tstep * cmat(1,1:ad%n_Z)
A_d = 1.d0 - theta * tstep * cmat(2,0:ad%n_Z)
A_u =      - theta * tstep * cmat(3,0:ad%n_Z-1)

! We have to do this manually because cmat is not a regular matrix
do i=1,ad%n_Z-1
  b(i)    = p(i) + (1.d0 - theta) * tstep * (p(i-1)*cmat(1,i) + p(i)*cmat(2,i) + p(i+1)*cmat(3,i))
enddo
b(0)      = p(0) + (1.d0 - theta) * tstep * (p(0)*cmat(2,i) + p(1)*cmat(3,i))
b(ad%n_Z) = p(ad%n_Z) + (1.d0 - theta) * tstep * (p(ad%n_Z-1)*cmat(1,ad%n_Z) + p(ad%n_Z)*cmat(2,ad%n_Z))

! Solve Ax=b for tridiagonal matrices. Result stored in b
call dgtsv(ad%n_Z+1,1,A_l,A_d,A_u,b,ad%n_Z+1,info)
if (info .ne. 0) write(*,*) 'info : ',info

! Output result
p = b
end subroutine coronal_timestep


!> Linear interpolation of coronal model charge at specific density and temperature
pure subroutine interpolate_coronal(cor, density, temperature, p_out, z_eff, rad)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_out !< distribution of charge states (sum = 1)
real*8, intent(out), optional   :: z_eff !< effective charge according to coronal equilibrium
real*8, intent(out), optional   :: rad !< radiated power according to coronal equilibrium

real*8, dimension(0:cor%n_Z)    :: Z !< The charge number at each charge state
integer                         :: iz

p = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature)

do iz=0,cor%n_Z
  if (p(iz)<0.) p(iz)=0.
end do

if (present(p_out)) then
  p_out = p
endif

if (present(z_eff)) then
  do iz=0,cor%n_Z
    Z(iz) = real(iz,8)
  enddo
  z_eff =  dot_product(p/sum(p),Z)
endif

if (present(rad)) then
  rad = L2Dinterp(cor%density,cor%temperature,cor%Prad(:,:),density,temperature)
endif
end subroutine interpolate_coronal


!> Linear interpolation of coronal model charge at specific density and temperature.
!> Evaluate the gradients only
pure subroutine interpolate_coronal_gradients(cor, density, temperature, p_Te_out, p_Ne_out, z_eff_Te, z_eff_Ne)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_Te_out, p_Ne_out !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, intent(out), optional   :: z_eff_Te, z_eff_Ne !< effective charge gradient according to coronal equilibrium

real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, dimension(0:cor%n_Z)    :: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, dimension(0:cor%n_Z)    :: Z !< The charge number at each charge state
integer                         :: iz

p    = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature)

p_Te = L2D2interp_grad(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature,1)
p_Ne = L2D2interp_grad(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature,2)

! Converting log gradient to real gradient
p_Te = p_Te / (log(10.)*10.0**temperature)
p_Ne = p_Ne / (log(10.)*10.0**density)

if (present(p_Te_out)) then
  p_Te_out = p_Te/sum(p)
endif

if (present(p_Ne_out)) then
  p_Ne_out = p_Ne/sum(p)
endif

if (present(z_eff_Te)) then
  do iz=0,cor%n_Z
    Z(iz) = real(iz,8)
  enddo
  z_eff_Te =  dot_product(p_Te/sum(p),Z)
endif

if (present(z_eff_Ne)) then
  do iz=0,cor%n_Z
    Z(iz) = real(iz,8)
  enddo
  z_eff_Ne =  dot_product(p_Ne/sum(p),Z)
endif
end subroutine interpolate_coronal_gradients

!> This is to output a coronal equilibrium charge distribution as a
!> function of temperature assuming constant density 10^20/m^3
!> to a file charge_distribution.dat
subroutine output_coronal(cor)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type

! Temporary variable for charge state distribution
integer             :: i_T, i_ion
real*8, allocatable :: dP_imp_dT(:), P_imp(:)
real*8              :: T_rad
real*8              :: Z_imp


open(20,file="charge_distribution.dat")

write(20,'(A11)',advance='no') 'temperature (log10(K))', 'charge states'
write(20,'(A11)') 'summation of all states'

do i_T = 0, 100

  T_rad = 2. + 0.06*real(i_T,8)

  if (allocated(P_imp)) deallocate(P_imp)
  if (allocated(dP_imp_dT)) deallocate(dP_imp_dT)

  allocate(P_imp(0:cor%n_Z))
  allocate(dP_imp_dT(0:cor%n_Z))
  call cor%interp(density=20.d0,temperature=T_rad,p_out=P_imp,z_eff=Z_imp)
  write(20,'(f12.3)',advance='no'), T_rad
  do i_ion = 0, cor%n_Z
    write(20,'(f12.5)',advance='no'), P_imp(i_ion)
  end do
  write(20,'(f12.5)'), sum(P_imp)
end do
close (20)
end subroutine output_coronal

end module mod_coronal

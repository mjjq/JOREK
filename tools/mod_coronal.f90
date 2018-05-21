!> module takes the OPEN-ADAS data to calculate the coronal equil-
!> ibrium temperature and radiation.
!> If you need time-dependent solutions of the corona matrix timestepping look
!> in the revision history for this file.
module mod_coronal
use mod_openadas
implicit none
private
public coronal
public output_coronal

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

real*8, dimension(0:ad%n_Z) :: p
integer :: n_d, n_T, iz, k, m
real*8 :: ion_rate, rec_rate

cor%n_Z = ad%n_Z
n_d = 10
n_T = 1000

allocate(cor%density(n_d), cor%temperature(n_T), cor%Z(n_d,n_T,0:cor%n_Z), cor%Prad(n_d,n_T))
do m=1, n_d
  cor%density(m) = 18.d0 + float(m-1)/n_d * (21.-18.) ! log10 [m^-3], linear between 18 and 21

  do k=1, n_T
    cor%temperature(k) = log10( 1.d0 + exp(log(4.d4)*float(k-1)/(float(n_T-1))) - 1.d0 ) + log10(EL_CHG) - log10(K_BOLTZ) ! in log10 [K]
    ! 1 to 40000 eV in logscale

    p(0) = 1.d0
    do iz=1,ad%n_Z
      ion_rate = ad%SCD%interp(iz, cor%density(m), cor%temperature(k)) ! ionizing to level iz (0 is neutral)
      rec_rate = ad%ACD%interp(iz, cor%density(m), cor%temperature(k)) ! recombining from iz+1
      p(iz) = p(iz-1) * ion_rate/rec_rate
    end do

    cor%Z(m,k,:)  = p/sum(p)
    cor%Prad(m,k) = coronal_Prad(ad, cor%density(m), cor%temperature(k), p/sum(p)) ! Do not set neutral density yet
  enddo
enddo
end function coronal_equilibrium


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
!> plot with gnuplot like
!> 
!> set logscale y
!> p for [i=2:20] 'charge_distribution.dat' u 1:i t ''.(i-2) w l
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

do i_T = 1, size(cor%temperature,1)
  T_rad = cor%temperature(i_T)

  if (allocated(P_imp)) deallocate(P_imp)
  if (allocated(dP_imp_dT)) deallocate(dP_imp_dT)

  allocate(P_imp(0:cor%n_Z))
  allocate(dP_imp_dT(0:cor%n_Z))
  call cor%interp(density=20.d0,temperature=T_rad,p_out=P_imp,z_eff=Z_imp)
  write(20,'(f12.3)',advance='no') T_rad
  do i_ion = 0, cor%n_Z
    write(20,'(f12.5)',advance='no') P_imp(i_ion)
  end do
  write(20,'(f12.5)') sum(P_imp)
end do
close (20)
end subroutine output_coronal

end module mod_coronal

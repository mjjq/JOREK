!> module takes the OPEN-ADAS data to calculate the coronal equil-
!> ibrium temperature and radiation.
!> If you need time-dependent solutions of the corona matrix timestepping look
!> in the revision history for this file.
module mod_coronal
use mod_openadas
use mod_interp_splinear
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
  type(Fspline)       :: ZFspline  !< Spline functions for effective charge
  type(Fspline)       :: PradFspline  !< Spline functions for CE radiation function
  type(Fspline), allocatable :: PFspline(:)  !< Spline functions for each charge state
contains
  procedure :: interp => interpolate_coronal
  procedure :: interp_gradients => interpolate_coronal_gradients
  procedure :: interp_spl => interpolate_coronal_spl
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
real*8, allocatable :: Z_eff(:,:)

cor%n_Z = ad%n_Z
n_d = 10
n_T = 1000

allocate(cor%density(n_d), cor%temperature(n_T), cor%Z(n_d,n_T,0:cor%n_Z), cor%Prad(n_d,n_T))
allocate(Z_eff(n_d,n_T))
Z_eff = 0.0

allocate(cor%PFspline(0:cor%n_Z))
do iz=0,ad%n_Z
  call AllocFspline(cor%PFspline(iz),n_T,n_d)
end do
call AllocFspline(cor%ZFspline,n_T,n_d)
call AllocFspline(cor%PradFspline,n_T,n_d)

do m=1, n_d
  cor%density(m) = 18.d0 + float(m-1)/n_d * (21.-18.) ! log10 [m^-3], linear between 18 and 21
end do
do k=1, n_T
  cor%temperature(k) = log10( 1.d0 + exp(log(4.d4)*float(k-1)/(float(n_T-1))) - 1.d0 ) + log10(EL_CHG) - log10(K_BOLTZ) ! in log10 [K], 1 to 40000 eV in logscale
end do

cor%ZFspline%xspline = cor%temperature
cor%ZFspline%ylinear = cor%density
cor%PradFspline%xspline = cor%temperature
cor%PradFspline%ylinear = cor%density
do iz=0,ad%n_Z
  cor%PFspline(iz)%xspline = cor%temperature
  cor%PFspline(iz)%ylinear = cor%density
end do

do m=1, n_d
  do k=1, n_T
    p(0) = 1.d0
    do iz=1,ad%n_Z
      ion_rate = ad%SCD%interp(iz, cor%density(m), cor%temperature(k)) ! ionizing to level iz (0 is neutral)
      rec_rate = ad%ACD%interp(iz, cor%density(m), cor%temperature(k)) ! recombining from iz+1
      p(iz) = p(iz-1) * ion_rate/rec_rate
    end do

    cor%Z(m,k,:)  = p/sum(p)
    do iz=1,ad%n_Z
      Z_eff(m,k) = Z_eff(m,k) + cor%Z(m,k,iz) * real(iz,8)
    end do
    cor%Prad(m,k) = coronal_Prad(ad, cor%density(m), cor%temperature(k), p/sum(p)) ! Do not set neutral density yet
  enddo

  cor%ZFspline%xspline = cor%temperature

  call spline(n_T,cor%ZFspline%xspline,Z_eff(m,:),0.d0,0.d0,2,&
                  cor%ZFspline%Aspline(m,:),cor%ZFspline%Bspline(m,:),&
                  cor%ZFspline%Cspline(m,:),cor%ZFspline%Dspline(m,:))

  call spline(n_T,cor%PradFspline%xspline,cor%Prad(m,:),0.d0,0.d0,2,&
                  cor%PradFspline%Aspline(m,:),cor%PradFspline%Bspline(m,:),&
                  cor%PradFspline%Cspline(m,:),cor%PradFspline%Dspline(m,:))
  do iz=0,ad%n_Z
    call spline(n_T,cor%PFspline(iz)%xspline,cor%Z(m,:,iz),0.d0,0.d0,2,&
                    cor%PFspline(iz)%Aspline(m,:),cor%PFspline(iz)%Bspline(m,:),&
                    cor%PFspline(iz)%Cspline(m,:),cor%PFspline(iz)%Dspline(m,:))
  end do
enddo

end function coronal_equilibrium


!> Linear interpolation of coronal model charge at specific density and temperature
subroutine interpolate_coronal(cor, density, temperature, p_out, z_eff, rad)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_out !< distribution of charge states (sum = 1)
real*8, intent(out), optional   :: z_eff !< effective charge according to coronal equilibrium
real*8, intent(out), optional   :: rad !< radiated power according to coronal equilibrium

real*8, dimension(0:cor%n_Z)    :: Z !< The charge number at each charge state
integer                         :: iz

!p = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature)
do iz = 0, cor%n_z
  call SL2Dinterp(cor%PFspline(iz),temperature,density,fout=p(iz))
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
  !rad = L2Dinterp(cor%density,cor%temperature,cor%Prad(:,:),density,temperature)
  call SL2Dinterp(cor%PradFspline,temperature,density,fout=rad)
endif
end subroutine interpolate_coronal


!> Linear interpolation of coronal model charge at specific density and temperature.
!> Evaluate the gradients only
subroutine interpolate_coronal_gradients(cor, density, temperature, p_Te_out, p_Ne_out, &
                                         z_eff_Te, z_eff_Ne, rad_Te, rad_Ne)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_Te_out, p_Ne_out !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, intent(out), optional   :: z_eff_Te, z_eff_Ne !< effective charge gradient according to coronal equilibrium
real*8, intent(out), optional   :: rad_Te, rad_Ne !< radiation function gradient

real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, dimension(0:cor%n_Z)    :: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, dimension(0:cor%n_Z)    :: Z !< The charge number at each charge state
integer                         :: iz

!p    = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature)
!p_Te = L2D2interp_grad(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature,1)
!p_Ne = L2D2interp_grad(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature,2)

do iz = 0, cor%n_z
  call SL2Dinterp(cor%PFspline(iz),temperature,density,fout=p(iz),dfout_dx=p_Te(iz),dfout_dy=p_Ne(iz))
  if (p(iz)<0.) p(iz)=0.
end do

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

if (present(rad_Te)) then
  call SL2Dinterp(cor%PradFspline,temperature,density,dfout_dx=rad_Te)
end if
if (present(rad_Ne)) then
  call SL2Dinterp(cor%PradFspline,temperature,density,dfout_dy=rad_Ne)
end if


end subroutine interpolate_coronal_gradients

! Spline-linear interpolation of the coronal equilibrium for both value and
! gradients
subroutine interpolate_coronal_spl(cor, density, temperature, p_out, p_Te_out, p_Ne_out, z_out,&
                                         z_Te_out, z_TeTe_out, z_Ne_out, rad_out, rad_Te_out, rad_Ne_out)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_out, p_Te_out, p_Ne_out !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, intent(out), optional   :: z_out, z_Te_out, z_TeTe_out, z_Ne_out !< effective charge gradient according to coronal equilibrium
real*8, intent(out), optional   :: rad_out, rad_Te_out, rad_Ne_out !< density multiplied radiation function and gradients

real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, dimension(0:cor%n_Z)    :: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
integer                         :: iz

real*8                          :: z, z_Te, z_Ne, z_TeTe ! local variables preparing for output
real*8                          :: rad, rad_Te, rad_Ne ! local variables preparing for output

if (present(p_out) .or. present(p_Te_out) .or. present(p_Ne_out)) then
  do iz = 0, cor%n_z
    call SL2Dinterp(cor%PFspline(iz),temperature,density,fout=p(iz),dfout_dx=p_Te(iz),dfout_dy=p_Ne(iz))
    if (p(iz)<0.) p(iz)=0.
  end do

  ! Converting log gradient to real gradient
  p_Te = p_Te / (log(10.)*10.0**temperature)
  p_Ne = p_Ne / (log(10.)*10.0**density)

  if (present(p_out)) p_out = p
  if (present(p_Te_out)) p_Te_out = p_Te/sum(p)
  if (present(p_Ne_out)) p_Ne_out = p_Ne/sum(p)
end if

if (present(z_out) .or. present(z_Te_out) .or. present(z_TeTe_out) .or. present(p_Ne_out)) then
  call SL2Dinterp(cor%ZFspline,temperature,density,fout=z,dfout_dx=z_Te,dfout_dy=z_Ne,d2fout_dx2=z_TeTe)

  ! Converting log gradient to real gradient
  z_Te = z_Te / (log(10.)*10.0**temperature)
  z_Ne = z_Ne / (log(10.)*10.0**density)
  z_TeTe = (z_TeTe / (log(10.)**2.0 * 10.0**(2.*temperature))) - z_Te/(10.0**temperature)

  if (present(z_out)) z_out = z
  if (present(z_Te_out)) z_Te_out = z_Te
  if (present(z_Ne_out)) z_Ne_out = z_Ne
  if (present(z_TeTe_out)) z_TeTe_out = z_TeTe
end if

if (present(rad_out) .or. present(rad_Te_out) .or. present(rad_Ne_out)) then
  call SL2Dinterp(cor%PradFspline,temperature,density,fout=rad,dfout_dx=rad_Te,dfout_dy=rad_Ne)

  ! Converting log gradient to real gradient
  rad_Te = rad_Te / (log(10.)*10.0**temperature)
  rad_Ne = rad_Ne / (log(10.)*10.0**density)

  if (present(rad_out)) rad_out = rad
  if (present(rad_Te_out)) rad_Te_out = rad_Te
  if (present(rad_Ne_out)) rad_Ne_out = rad_Ne
end if

end subroutine interpolate_coronal_spl

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
real*8, allocatable :: P_imp(:)
real*8              :: T_rad, Lrad
real*8              :: Z_eff

open(20,file="charge_distribution.dat")

write(20,'(4A22)',advance='no') 'temperature (log10(K))', 'charge states', 'summation', 'effective charge'
write(20,'(A22)') 'radiation function'

do i_T = 1, size(cor%temperature,1)
  T_rad = cor%temperature(i_T)

  if (allocated(P_imp)) deallocate(P_imp)

  allocate(P_imp(0:cor%n_Z))
  call cor%interp_spl(density=20.d0,temperature=T_rad,p_out=P_imp,z_out=Z_eff,rad_out=Lrad)
  Lrad = Lrad / (1.d20) ! This is to recover the radiation coefficient
  write(20,'(f12.3)',advance='no') T_rad
  do i_ion = 0, cor%n_Z
    write(20,'(f12.5)',advance='no') P_imp(i_ion)
  end do
  write(20,'(f12.5)',advance='no') sum(P_imp)
  write(20,'(f12.5)',advance='no') Z_eff
  write(20,'(e14.6)') Lrad
end do
close (20)

end subroutine output_coronal

end module mod_coronal

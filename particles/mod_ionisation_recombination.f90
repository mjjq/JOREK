!> This module contains routines to calculate ionisation
!> and recombination probabilities of particles in a specific time
module mod_ionisation_recombination
use mod_openadas
implicit none
private
public new_charge
public fields_interp_ne_Te

contains

!> Calculate new charge state at a specific density, temperature and timestep
pure function new_charge(z, ad, electron_density, electron_temperature, timestep, ran2) result(z_new)
implicit none

integer, intent(in)              :: z !< Old charge state
type(ADF11_all), intent(in)      :: ad !< ADF11 data
real*8, intent(in)               :: electron_density !< log10 Electron density in m^-3
real*8, intent(in)               :: electron_temperature !< log10 Electron temperature in K
real*8, intent(in)               :: timestep !< Timestep in s
real*8, intent(in), dimension(2) :: ran2 !< Random numbers to use to select a new charge or not
integer :: z_new

real*8 :: prob(2)

z_new = z
! probabilities of recombination and ionisation events
prob = 1.d0 - exp(-[ad%ACD%interp(z,   electron_density, electron_temperature), & ! rec
                    ad%SCD%interp(z+1, electron_density, electron_temperature)] & ! ion
                  * 10.d0**electron_density * timestep)
z_new = z_new + sum([-1, 1], mask=(prob .gt. ran2))

end function new_charge

pure subroutine fields_interp_ne_Te(fields, time, s, t, phi, i_elm, n_e, T_e)
use mod_fields
use phys_module, only: central_density
use constants
class(fields_base), intent(in)                    :: fields
real*8, intent(in)                                :: time, s, t, phi
integer, intent(in)                               :: i_elm
real*8, intent(out)                               :: n_e !< electron density [m^-3]
real*8, intent(out)                               :: T_e !< electron temperature [eV]

real*8, dimension(2) :: P, P_s, P_t, P_phi, P_time
real*8               :: R, R_s, R_t, Z, Z_s, Z_t
call fields%interp_PRZ(time,i_elm,&
#if (JOREK_MODEL == 400)
      [5,8],& ! electron temperature
#else
      [5,6],& ! electron temperature + ion temperature (assumed equal)
#endif
          2,s,t,phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)

n_e = P(1) * 1d20                           ! plasma density [1/m^3]
T_e = P(2)/(2.d0*MU_ZERO*central_density*1.d20)/K_BOLTZ
#if (JOREK_MODEL == 400)
T_e = T_e*2d0 ! P(1) contains the electron temperature, reverse previous correction
#endif
end subroutine fields_interp_ne_Te
end module mod_ionisation_recombination

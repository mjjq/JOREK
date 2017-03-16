!> This module contains routines to calculate ionisation
!> and recombination probabilities of particles in a specific time
module mod_ionisation_recombination
use mod_openadas
implicit none
private
public new_charge

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
end module mod_ionisation_recombination

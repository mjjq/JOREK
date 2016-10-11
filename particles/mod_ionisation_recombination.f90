!> This module contains routines to calculate ionisation
!> and recombination probabilities of particles in a specific time
module mod_ionisation_recombination
use mod_openadas
implicit none
private
public new_charge

contains

!> Calculate new charge state at a specific density, temperature and timestep
!> TODO remove dependency on random_number and create a pure function
function new_charge(z, ad, electron_density, electron_temperature, timestep) result(z_new)
implicit none

integer, intent(in)              :: z !< Old charge state
type(type_ADF11_all), intent(in) :: ad !< ADF11 data
real*8, intent(in)               :: electron_density !< log10 Electron density in m^-3
real*8, intent(in)               :: electron_temperature !< log10 Electron temperature in K
real*8, intent(in)               :: timestep !< Timestep in s
integer :: z_new

real*8 :: prob(2)
real*8 :: rand(2)

z_new = z
call random_number(rand)
! probabilities of recombination and ionisation events
prob = 1.d0 - exp([GRC(ad%ACD, z,   electron_density, electron_temperature), & ! rec
                   GRC(ad%SCD, z+1, electron_density, electron_temperature)] & ! ion
         * 10.d0**electron_density * timestep)
if (prob(1) .gt. rand(1)) z_new = z_new - 1 ! recombination
if (prob(2) .gt. rand(2)) z_new = z_new + 1 ! ionization

end function new_charge
end module mod_ionisation_recombination

!> This module contains routines to calculate ionisation
!! and recombination probabilities of particles in a specific time
module mod_ionisation_recombination
use openadas

contains
function new_charge(z, ad, electron_density, electron_temperature, timestep) result(z_new)
implicit none

integer, intent(in)              :: z !< Old charge state
type(type_ADF11_all), intent(in) :: ad !< ADF11 data
real*8, intent(in)               :: electron_density
real*8, intent(in)               :: electron_temperature
real*8, intent(in)               :: timestep
integer :: z_new

real*8 :: prob(2)
real*8 :: rand(2)

z_new = z
call random_number(rand)
! probabilities of recombination and ionisation events
prob = (/GRC(ad%ACD, z, electron_density, electron_temperature), & ! rec
         GRC(ad%SCD, z+1,   electron_density, electron_temperature)/) * 10.d0**electron_density * timestep ! ion
if (prob(1) .gt. rand(1)) z_new = z_new - 1 ! recombination
if (prob(2) .gt. rand(2)) z_new = z_new + 1 ! ionization
if (prob(1) .gt. 1.d0 .and. prob(2) .gt. 1.d0) write(*,*) "WARNING: ion AND rec probabilities >1, lower timestep!", prob

end function new_charge
end module mod_ionisation_recombination

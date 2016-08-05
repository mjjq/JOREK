!> Module containing mathematical and physical constants which are used in the code
module mod_constants
  implicit none

  real*8,  parameter :: PI               = 3.1415926535897932385d0  !< Half-circle constant
  real*8,  parameter :: TWOPI            = 6.2831853071795864769d0  !< Circle constant
  real*8,  parameter :: MU_ZERO          = 4.d-7*PI                 !< Magnetic constant  [Vs/Am]
  real*8,  parameter :: EPS_ZERO         = 8.854187817d-12          !< Vacuum permittivity [F/m]
  real*8,  parameter :: EL_CHG           = 1.602176565d-19          !< Elementary charge  [C]
  real*8,  parameter :: K_BOLTZ          = 1.3806488d-23            !< Boltzmann constant [J/K]
  real*8,  parameter :: MASS_PROTON      = 1.67262178d-27           !< proton mass [kg]
  real*8,  parameter :: ATOMIC_MASS_UNIT = 1.660538921d-27          !< standardized mass unit [kg]
  real*8,  parameter :: MASS_ELECTRON    = 9.10938291d-31           !< electron mass [kg]
end module mod_constants

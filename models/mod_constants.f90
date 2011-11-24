!> Module containing constants which are used in the code
module constants
  
  implicit none
  
  ! @name Mathematical and physical constants
  real*8,  parameter :: PI         = 3.1415926535897932385d0
  real*8,  parameter :: MU_ZERO    = 4.d-7*PI                 !< Magnetic constant  [Vs/Am]
  real*8,  parameter :: EL_CHG     = 1.602176565d-19          !< Elementary charge  [C]
  real*8,  parameter :: K_BOLTZ    = 1.3806488d-23            !< Boltzmann constant [J/K]
  
  !> @name Constants which describe the domain of a certain position (used by function which_domain)
  integer, parameter :: DOMAIN_PLASMA         = 0    !< Plasma region
  integer, parameter :: DOMAIN_SOL            = 1    !< Scrape-off layer
  integer, parameter :: DOMAIN_OUTER_SOL      = 2    !< Outer scrape-off layer (double-null)
  integer, parameter :: DOMAIN_UPPER_PRIVATE  = 3    !< Upper private flux region
  integer, parameter :: DOMAIN_LOWER_PRIVATE  = 4    !< Lower private flux region
  
  !> @name Parameters which describe the X-point case
  integer, parameter :: LOWER_XPOINT          = 1
  integer, parameter :: UPPER_XPOINT          = 2
  integer, parameter :: DOUBLE_NULL           = 3
  
end module constants

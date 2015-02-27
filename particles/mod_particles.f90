module particles

use parameters

type type_particle

  real*8  :: x(3)             !< particle position in real space (R,Z,phi)
  real*8  :: st(n_dim)        !< particle position in the finite element (i_elm)
  real*8  :: v(3)             !< particle velocity in (R,Z,phi)
  real*8  :: q
  real*8  :: mass
  real*8  :: weight
  real*8  :: temperature
  real*8  :: radiation
  integer :: i_elm

end type type_particle
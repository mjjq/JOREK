!> This is essentially the same as mod_integrals3D, but without MPI. It is meant
!! for use in diagnostic programs like jorek2_postproc.
module mod_integrals3D_nompi

#define NOMPIVERSION=1

#include "mod_integrals3D.f90"

#undef NOMPIVERSION

end module mod_integrals3D_nompi

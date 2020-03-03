!> This is essentially the same as mod_integrals3D, but without MPI. It is meant
!! for use in diagnostic programs like jorek2_postproc.
module mod_integrals3D_nompi

use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use tr_module
use phys_module
use pellet_module
use mod_interp
use convert_character
use mod_expression
use mod_resistivity
use mod_poloidal_currents, only : integrated_normal_bnd_curr 
use corr_neg
use equil_info, only : get_psi_n, ES

#define NOMPIVERSION=1
#include "mod_integrals3D.f90"
#undef NOMPIVERSION


end module mod_integrals3D_nompi

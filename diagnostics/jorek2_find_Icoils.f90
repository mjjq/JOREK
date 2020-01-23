!> Program to find PF coil currents given a coil geometry file (coils_geo.txt) 
!> The found Icoils try to match the external field of the jorek_restart.h5 plasma
program jorek2_find_Icoils

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use mpi_mod
use mod_import_restart
use mod_vtk
use mod_interp
use mod_plasma_response
use vacuum
use mod_boundary

implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list), pointer :: bnd_node_list    

integer                   :: i, my_id, ifail, ierr
character                 :: buffer*80, bnd_file*20

type(t_coil), allocatable :: coils(:)
integer                   :: n_pt 
real*8, allocatable       :: R(:), Z(:), psi_p(:), psi_c(:)


call flush_it(6)

allocate(node_list)
allocate(element_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

call flush_it(6)

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)

call initialise_basis                              ! define the basis functions at the Gaussian points
call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

! --- Read coils geometry given by filaments
call read_coils(coils, 'coils_geo.txt')

! --- Find PF coils currents
if (tokamak_device == 'JET') then
  call find_Icoils_JET(node_list,element_list,bnd_node_list,bnd_elm_list, coils)
else
  call find_Icoils(node_list,element_list,bnd_node_list,bnd_elm_list, coils)
endif

end program jorek2_find_Icoils

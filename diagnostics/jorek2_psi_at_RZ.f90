!> Program to calculate psi in an arbitrary set of (R,Z) points
!> Useful to calculate psi outside the JOREK's domain
!> It needs a coil geometry file with the proper PF coil currents in
!> the input file (pf_coils(:)%current)
!> If you don't know the coil currents use jorek2_find_Icoils
program jorek2_psi_at_RZ

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

bnd_file="RZ_bnd.dat"

! --- Read R,Z points in which the total flux will be calculated 
open(42, file=trim(bnd_file), action='read', status='old', iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'ERROR: Cannot open RZ file "'//trim(bnd_file)//'".'
    stop
  end if
  read(42,*) n_pt
  allocate(R(n_pt), Z(n_pt), psi_p(n_pt), psi_c(n_pt))
  do i=1, n_pt
    read(42,*) R(i), Z(i)
  enddo
close(42)

! --- Read coils geometry given by filaments
call read_coils(coils, 'coils_geo.txt')

! --- Calculate flux produced by the plasma current
call psi_plasma(node_list,element_list, R, Z, psi_p)

! --- Calculate the flux produced by the coil currents
call psi_coils (coils, R, Z, psi_c)

! --- Output psi in the jorek input file format
open(43, file='RZ_psi_input.dat', action='write', iostat=ierr)
  do i=1, n_pt
    write(43,'(A,i3,A,e16.8,A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
             '  R_boundary(',i,') =',R(i), &
             ', Z_boundary(',i,') =',Z(i), &
             ', psi_boundary(',i,') =',psi_p(i)+psi_c(i),','
  enddo
close(43)


end program jorek2_psi_at_RZ

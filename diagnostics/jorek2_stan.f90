!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2vtk

use parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use nodes_elements

implicit none

type (type_surface_list) 		:: flux_list
integer					:: i, k_tor, my_id, ierr, ifail
integer					:: i_elm_xpoint(2), i_elm_axis
real*8					:: R_xpoint(2), Z_xpoint(2), psi_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8					:: R_axis,      Z_axis,      psi_axis,      s_axis,      t_axis
real*8					:: psi_bnd

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list,element_list, 'jorek_restart.rst', rst_format, ierr)

call initialise_basis                              ! define the basis functions at the Gaussian points

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

if (xpoint) then
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
  psi_bnd  = psi_xpoint(1)
  if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
    psi_bnd = psi_xpoint(2)
  endif
else
  psi_bnd = 0.d0
endif




flux_list%n_psi = 100
allocate(flux_list%psi_values(flux_list%n_psi))
do i=1,flux_list%n_psi-1
  flux_list%psi_values(i) = psi_axis + 1.2 * (psi_bnd - psi_axis) * real(i) / real(flux_list%n_psi-1)
enddo

flux_list%psi_values(flux_list%n_psi) = psi_bnd

call find_flux_surfaces(xpoint,xcase,node_list,element_list,flux_list)

call reorder_flux_surfaces(node_list, element_list, flux_list, ierr)
write(*,*)'Error from reorder_flux_surfaces:',ierr



deallocate(flux_list%psi_values)
return

 
end program jorek2vtk

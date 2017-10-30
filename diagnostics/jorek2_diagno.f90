!**********************************************************************
!* program to extract data from a JOREK2 restart file                 *
!**********************************************************************

program jorek2_diagno
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
use pellet_module
use mpi_mod
use mod_import_restart
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
integer :: i, in, i_tor
real*8  :: growth_kin, growth_mag,density,density_in,density_out,pressure,pressure_in,pressure_out
real*8  :: Rplot(2), Zplot(2)
real*8  :: psi_axis,R_axis,Z_axis,s_axis,t_axis
integer :: ifail, my_id, ierr, i_elm_axis
integer :: required, provided, StatInfo


write(*,*) '***************************************'
write(*,*) '* JOREK2_diagno                       *'
write(*,*) '***************************************'

my_id=0

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
call MPI_Init_thread(required, provided, StatInfo)


call initialise_parameters(my_id, "__NO_FILENAME__")

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr)

call initialise_basis                              ! define the basis functions at the Gaussian points

call Integrals_3D(my_id,node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

end program jorek2_diagno

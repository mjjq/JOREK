!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2_get_new_extension_patch_connections

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp
use pellet_module
use mpi_mod
use mod_import_restart

implicit none

type (type_node_list)   , pointer :: node_list
type (type_element_list), pointer :: element_list

integer :: my_id, k_tor, i_node, ierr, i_ext, i_part
integer :: i_bnd_beg, i_bnd_end
real*8  :: diff_min_beg, diff_min_end, diff


write(*,*) '***************************************'
write(*,*) '*       jorek2vtk                     *'
write(*,*) '***************************************'

allocate(node_list)
allocate(element_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)
call initialise_basis                              ! define the basis functions at the Gaussian points

! --- First, find out which bnd nodes are our starting/ending points
do i_ext = 1,n_wall_blocks
  diff_min_beg = 1.d10
  diff_min_end = 1.d10
  do i_node = 1,node_list%n_nodes
    if (node_list%node(i_node)%boundary .eq. 0) cycle
    diff = sqrt( (node_list%node(i_node)%x(1,1)-R_block_points_left(i_ext,1))**2 &
                +(node_list%node(i_node)%x(1,2)-Z_block_points_left(i_ext,1))**2 )
    if (diff .lt. diff_min_beg) then
      diff_min_beg = diff
      i_bnd_beg = i_node
    endif
    diff = sqrt( (node_list%node(i_node)%x(1,1)-R_block_points_right(i_ext,1))**2 &
                +(node_list%node(i_node)%x(1,2)-Z_block_points_right(i_ext,1))**2 )
    if (diff .lt. diff_min_end) then
      diff_min_end = diff
      i_bnd_end = i_node
    endif
  enddo
  R_block_points_left (i_ext,1) = node_list%node(i_bnd_beg)%x(1,1)
  Z_block_points_left (i_ext,1) = node_list%node(i_bnd_beg)%x(1,2)
  R_block_points_right(i_ext,1) = node_list%node(i_bnd_end)%x(1,1)
  Z_block_points_right(i_ext,1) = node_list%node(i_bnd_end)%x(1,2)
  write(*,'(A,i2,A,i2)')            ' n_ext_block         (',i_ext,           ')    = ',n_ext_block(i_ext)
  write(*,'(A,i2,A,i2)')            ' n_block_points_left (',i_ext,           ')    = ',n_block_points_left(i_ext)
  do i_part=1,n_block_points_left(i_ext)
    write(*,'(A,i2,A,i1,A,sp,f7.4)')' R_block_points_left (',i_ext,',',i_part,')  = ',R_block_points_left (i_ext,i_part)
    write(*,'(A,i2,A,i1,A,sp,f7.4)')' Z_block_points_left (',i_ext,',',i_part,')  = ',Z_block_points_left (i_ext,i_part)
  enddo
  write(*,'(A,i2,A,i1)')            ' n_block_points_right(',i_ext,           ')    = ',n_block_points_right(i_ext)
  do i_part=1,n_block_points_right(i_ext)
    write(*,'(A,i2,A,i1,A,sp,f7.4)')' R_block_points_right(',i_ext,',',i_part,')  = ',R_block_points_right(i_ext,i_part)
    write(*,'(A,i2,A,i1,A,sp,f7.4)')' Z_block_points_right(',i_ext,',',i_part,')  = ',Z_block_points_right(i_ext,i_part)
  enddo
  write(*,*)''
enddo

end program jorek2_get_new_extension_patch_connections








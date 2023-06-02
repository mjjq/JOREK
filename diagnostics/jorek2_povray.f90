!> 3D visualizations of JOREK data with the raytracing program Povray.
!!
!! 2023 by Matthias Hoelzl and Nina Schwarz
!!
program jorek2_povray

use constants
use data_structure
use nodes_elements
use phys_module
use mod_import_restart
use mod_interp
use basis_at_gaussian
use mod_basisfunctions
use mod_boundary
use equil_info

implicit none

integer, parameter :: ifile   = 53 !< File handle
logical, parameter :: surface = .true.
real*8,  parameter :: phimin  = 0.d0
real*8,  parameter :: phimax  = 270.d0 / 360.d0 * 2.d0 * PI
integer, parameter :: n_phi   = 64
integer, parameter :: nsub    = 4

integer :: ierr, my_id, i_elm, i_s, i_t, iv, idof
real*8  :: s, t, phi
real*8  :: HH(nsub+1,nsub+1,4,n_degrees)
real*8  :: val(nsub+1,nsub+1)
real*8  :: R(nsub+1,nsub+1)
real*8  :: Z(nsub+1,nsub+1)
real*8  :: x(nsub+1,nsub+1)
real*8  :: y(nsub+1,nsub+1)
type(type_element) :: element
type(type_node)    :: nodes(4)

! --- Initialize
my_id = 0
call initialise_parameters(my_id, "__NO_FILENAME__")
call det_modes()
call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
call initialise_basis()

call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)

call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)

open(ifile, file='jorek_3d.pov', form='formatted', status='replace')

do i_s = 1, nsub-1
  s = real(i_s-1)/real(nsub-1)
  do i_t = 1, nsub-1
    t = real(i_t-1)/real(nsub-1)
    call basisfunctions(s, t, HH(i_s, i_t, :, :))
  end do
end do

do i_elm = 1, element_list%n_elements
  element  = element_list%element(i_elm)
  nodes(:) = node_list%node(element%vertex(:))
  val(:,:) = 0.d0
  R(:,:)   = 0.d0
  Z(:,:)   = 0.d0
  do i_s = 1, nsub+1
    do i_t = 1, nsub+1
      do iv = 1, 4
        do idof = 1, 4
          R(i_s,i_t) = R(i_s,i_t) + nodes(iv)%x(1,idof,1)
          Z(i_s,i_t) = Z(i_s,i_t) + nodes(iv)%x(1,idof,2)
          !### val
        end do
      end do
    end do
  end do
  
  phi = 0.d0 !###
  x(:,:) = R(:,:) * sin(phi)
  y(:,:) = R(:,:) * cos(phi)
  
  
  do i_s = 1, nsub
    do i_t = 1, nsub
      call write_triangle_pov(ifile, &
        (/ x(i_s,i_t), x(i_s+1,i_t), x(i_s+1,i_t+1) /), &
        (/ y(i_s,i_t), y(i_s+1,i_t), y(i_s+1,i_t+1) /), &
        (/ z(i_s,i_t), z(i_s+1,i_t), z(i_s+1,i_t+1) /), &
        (/ 1.d0, 1.d0, 1.d0, 0.d0 /))
      call write_triangle_pov(ifile, &
        (/ x(i_s,i_t), x(i_s,i_t+1), x(i_s,i_t+1) /), &
        (/ y(i_s,i_t), y(i_s,i_t+1), y(i_s,i_t+1) /), &
        (/ z(i_s,i_t), z(i_s,i_t+1), z(i_s,i_t+1) /), &
        (/ 1.d0, 1.d0, 1.d0, 0.d0 /))
    end do
  end do
  
end do

close(ifile)



contains



subroutine write_triangle_pov(ifile, x, y, z, rgbt)
  integer,              intent(in) :: ifile
  real*8, dimension(3), intent(in) :: x, y, z
  real*8, dimension(4), intent(in) :: rgbt
  
  990 format(a)
  991 format('  <',f15.6,',',f15.6,',',f15.6,'>, <',f15.6,',',f15.6,',',f15.6,'>, <',f15.6, &
    ',',f15.6,',',f15.6,'>')
  992 format('  pigment{color rgbt<',f15.6,','f15.6,','f15.6,','f15.6,'>}')
  write(ifile,990) 'triangle {'
  write(ifile,991) x(1), y(1), z(1), x(2), y(2), z(2), x(3), y(3), z(3)
  write(ifile,992) rgbt
  write(ifile,990) '}'

end subroutine write_triangle_pov



end program jorek2_povray

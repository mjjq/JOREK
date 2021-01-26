!> Program to convert a JOREK2 restart file into binary VTK format
program compare_methods

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use mpi_mod
use mod_import_restart
use equil_info
use mod_boundary
use mod_interp
use mod_newton_methods

implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list 

integer :: my_id, i_elm, k_tor, ierr
real*8  :: psifind
real*8  :: psimin,  psimax
real*8  :: psimin2, psimax2
real*8  :: diff
real*8  :: Rmin,  Rmax
real*8  :: Zmin,  Zmax
real*8  :: Rmin2, Rmax2
real*8  :: Zmin2, Zmax2
real*8  :: diff_RZ
integer :: iv, im, is, n1, n2, ifail
real*8  :: p1, p2, p3, p4, dp1, dp4
real*8  :: a0, a1, a2, a3
real*8  :: st1, st2, st3
integer :: n_found, n_found2
real*8  :: st_found(n_order), st_found2(n_order)

allocate(node_list)
allocate(element_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                              ! define the basis functions at the Gaussian points

!call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

do i_elm=1,element_list%n_elements

  ! --- Check minmax functions
  if (.false.) then
    call psi_minmax(node_list,element_list,i_elm,psimin,psimax)
    call find_variable_minmax(node_list,element_list,i_elm, var_psi, psimin2, psimax2)
    diff = abs(psimin-psimin2) + abs(psimax-psimax2)
    call RZ_minmax(node_list,element_list,i_elm,Rmin,Rmax,Zmin,Zmax)
    call find_variable_minmax(node_list,element_list,i_elm, -1, Rmin2, Rmax2)
    call find_variable_minmax(node_list,element_list,i_elm, -2, Zmin2, Zmax2)
    diff_RZ = abs(Rmin-Rmin2) + abs(Rmax-Rmax2) + abs(Zmin-Zmin2) + abs(Zmax-Zmax2)
    if (diff    .gt. 1.d-15) write(*,'(A,i6,3e)')'psiminmax on element:',i_elm,diff,psimin,psimax
    if (diff_RZ .gt. 1.d-15) write(*,'(A,i6,1e)')'RZ minmax on element:',i_elm,diff_RZ
  endif

  ! --- Check root-solver
  if (.true.) then

    ! --- check psi_minmax
    call find_variable_minmax(node_list,element_list,i_elm, var_psi, psimin, psimax)
    psifind = 0.5 * (psimin+psimax)

    ! --- Do the 4 sides
    do iv=1, 4

      ! --- The old way, with cubic root formula
      im = MOD(iv,4) + 1
      n1 = element_list%element(i_elm)%vertex(iv)
      n2 = element_list%element(i_elm)%vertex(im)

      if (node_list%node(n1)%axis_node .and. node_list%node(n2)%axis_node) cycle

      is = mod(iv+1,2) + 2

      p1  =  node_list%node(n1)%values(1,1,1)  * element_list%element(i_elm)%size(iv,1)
      dp1 =  node_list%node(n1)%values(1,is,1) * element_list%element(i_elm)%size(iv,is)
      p4  =  node_list%node(n2)%values(1,1,1)  * element_list%element(i_elm)%size(im,1)
      dp4 =  node_list%node(n2)%values(1,is,1) * element_list%element(i_elm)%size(im,is)

      p2  = p1 + dp1
      p3  = p4 + dp4

      a3 = -        p1 + 3.d0 * p2 - 3.d0 * p3 + p4
      a2 = + 3.d0 * p1 - 6.d0 * p2 + 3.d0 * p3
      a1 = - 3.d0 * p1 + 3.d0 * p2
      a0 =          p1                              - psifind

      call SOLVP3(a0,a1,a2,a3,st1,st2,st3,ifail)

      n_found = 0
      if ((st1 .ge. 0.d0) .and. (st1 .le. 1.d0)) then
        n_found = n_found + 1
        st_found(n_found) = st1
      endif
      if ((st2 .ge. 0.d0) .and. (st2 .le. 1.d0)) then
        n_found = n_found + 1
        st_found(n_found) = st2
      endif
      if ((st3 .ge. 0.d0) .and. (st3 .le. 1.d0)) then
        n_found = n_found + 1
        st_found(n_found) = st3
      endif


      ! --- The new way, with Newton methods
      call newton_1D_find_value_on_element_side(node_list,element_list,i_elm, iv, var_psi, psifind, n_found2, st_found2)

      if (n_found .ne. n_found2) then
        write(*,'(A,4i6)') 'not found:',i_elm,iv,n_found,n_found2
      endif

      if (n_found2 .gt. 1) then
        write(*,'(A,3i6)') 'found many:',i_elm,iv,n_found2
      endif

      if ( (n_found .eq. 1) .and. (n_found2 .eq. 1) ) then
        ! --- IMPORTANT!!! find_flux_surfaces does it reversed for sides 3 and 4
        if (iv .ge. 3) st_found2(1) = 1.d0 - st_found2(1)
        ! --- Check difference
        diff = abs(st_found(1)-st_found2(1))
        if (diff .gt. 1.d-10) write(*,'(A,2i6,3e)')'found different',i_elm,iv,diff,st_found(1),st_found2(1)
        !write(*,'(A,2i6,3e)')'found???',i_elm,iv,diff,st_found(1),st_found2(1)
      endif

    enddo
  endif


enddo  ! n_elements


end program compare_methods

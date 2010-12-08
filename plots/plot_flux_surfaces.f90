subroutine plot_flux_surfaces(node_list,element_list,surface_list,frame,every_nth)
use data_structure
implicit none

! --- Routine parameters
type (type_node_list),    intent(in) :: node_list
type (type_element_list), intent(in) :: element_list
type (type_surface_list), intent(in) :: surface_list
logical,                  intent(in) :: frame
integer,                  intent(in) :: every_nth     ! Plot only every_nth flux surface

! --- internal variables
integer            :: i, j, k,ip, nplot, node1, node2, node3, node4, i_elm
real*8             :: t, rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi
real*8             :: R_min, R_max, Z_min, Z_max
real*8             :: dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9, dummy10
real*8,allocatable :: rplot(:), zplot(:)
character*13       :: LABEL

nplot=11
allocate(rplot(nplot),zplot(nplot))

LABEL= 'Flux surfaces'

R_min = 1.d20; R_max = -1.d20; Z_min  = 1.d20; Z_max = -1.d20
do i=1,node_list%n_nodes
  R_min = min(R_min,node_list%node(i)%x(1,1))
  R_max = max(R_max,node_list%node(i)%x(1,1))
  Z_min = min(Z_min,node_list%node(i)%x(1,2))
  Z_max = max(Z_max,node_list%node(i)%x(1,2))
enddo

if (frame) CALL NFRAME(21,11,1,R_min,R_max,Z_min,Z_max,LABEL,13,'R [m]',5,'Z [m]',5)

!call plot_grid(node_list,element_list,.true.,.false.)                               ! plot the grid

do j = 1, surface_list%n_psi, every_nth

!  write(*,*) ' plot : ',j,surface_list%flux_surfaces(j)%n_pieces

  do k=1,surface_list%flux_surfaces(j)%n_pieces

    i_elm = surface_list%flux_surfaces(j)%elm(k)

    node1 = element_list%element(i_elm)%vertex(1)
    node2 = element_list%element(i_elm)%vertex(2)
    node3 = element_list%element(i_elm)%vertex(3)
    node4 = element_list%element(i_elm)%vertex(4)

    rr1  = surface_list%flux_surfaces(j)%s(1,k)
    drr1 = surface_list%flux_surfaces(j)%s(2,k)
    rr2  = surface_list%flux_surfaces(j)%s(3,k)
    drr2 = surface_list%flux_surfaces(j)%s(4,k)

    ss1  = surface_list%flux_surfaces(j)%t(1,k)
    dss1 = surface_list%flux_surfaces(j)%t(2,k)
    ss2  = surface_list%flux_surfaces(j)%t(3,k)
    dss2 = surface_list%flux_surfaces(j)%t(4,k)

!    write(*,'(A,i5,4f10.4)') ' element : ',i_elm,rr1,ss1,rr2,ss2

    do ip = 1, nplot

      t = -1. + 2.*float(ip-1)/float(nplot-1)

      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

!      call INTERP2(RR(:,node1),RR(:,node2),RR(:,node3),RR(:,node4),ri,si,rplot(ip),dummy1,dummy2)
!      call INTERP2(ZZ(:,node1),ZZ(:,node2),ZZ(:,node3),ZZ(:,node4),ri,si,zplot(ip),dummy1,dummy2)

      call interp_RZ(node_list,element_list,i_elm,ri,si,rplot(ip),dummy1,dummy2,dummy3,dummy4,dummy5, &
                                                        zplot(ip),dummy6,dummy7,dummy8,dummy9,dummy10)

    enddo

    call lincol(1)
    if (surface_list%n_psi .eq. 1) call lincol(3)

    write(51,*) ' .5 setlinewidth '
    call lplot6(21,11,rplot,zplot,-nplot,' ')

  enddo

enddo

call lincol(0)

deallocate(rplot,zplot)

return
end subroutine plot_flux_surfaces

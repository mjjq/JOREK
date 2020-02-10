module py_plots_grids
  implicit none
contains







!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine py_plot_surface(filename,node_list,element_list,respline, ordered, surface_list, i_surf)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in) :: filename
  type (type_node_list),    intent(in) :: node_list
  type (type_element_list), intent(in) :: element_list
  logical,                  intent(in) :: respline
  logical,                  intent(in) :: ordered
  type (type_surface_list), intent(in) :: surface_list
  integer,                  intent(in) :: i_surf
  
  call print_py_plot_prepare_plot(filename)
  if (respline) then
    call respline_flux_surfaces(node_list,element_list,surface_list)
  endif
  if (ordered) then
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list, 'r', .false.)
  else
    call print_py_plot_unordered_flux_surfaces(filename, node_list, element_list, surface_list, i_surf)
  endif
  call print_py_plot_wall(filename)
  call print_py_plot_finish_plot(filename)


end subroutine py_plot_surface









!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine print_py_plot_prepare_plot(filename)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)          :: filename
  
  open(101,file=filename)
    write(101,'(A)')                          '#!/usr/bin/env python'
    write(101,'(A)')                          'import numpy as N'
    write(101,'(A)')                          'import pylab'
    write(101,'(A)')                          'def main():'
  close(101)



end subroutine print_py_plot_prepare_plot






  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine print_py_plot_finish_plot(filename)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)          :: filename
  
  open(101,file=filename,position='append')
    write(101,'(A)')                          ' pylab.axis("equal")'
    write(101,'(A)')                          ' pylab.show()'
    write(101,'(A)')                          ' '
    write(101,'(A)')                          'main()'
  close(101)



end subroutine print_py_plot_finish_plot






  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine print_py_plot_unordered_flux_surfaces(filename, node_list, element_list, surface_list, i_surf)

  use data_structure
  use mod_interp, only: interp_RZ
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in) :: filename
  type (type_node_list),    intent(in) :: node_list
  type (type_element_list), intent(in) :: element_list
  type (type_surface_list), intent(in) :: surface_list
  integer,                  intent(in) :: i_surf
  
  ! --- Internal variables
  integer       :: i, j, k, n_sub
  integer       :: i_elm
  integer       :: i_start, i_stop
  real*8        :: rr,    ss
  real*8        :: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8        :: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  real*8        :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2
  real*8        :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8        :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8        :: RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss
  real*8        :: ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss
  real*8        :: dRRg1_dt, dZZg1_dt, dRRg2_dt, dZZg2_dt
  
  n_sub = 2 ! 2 minimum!
  
  i_start = 1
  i_stop  = surface_list%n_psi
  if (i_surf .gt. 0) then
    i_start = i_surf
    i_stop  = i_surf
  endif
  
  open(101,file=filename,position='append')
    write(101,'(A,i2,A)')                     ' rplot = N.zeros(',n_sub,')'
    write(101,'(A,i2,A)')                     ' zplot = N.zeros(',n_sub,')'
    do i= i_start, i_stop
      do j=1,surface_list%flux_surfaces(i)%n_pieces
        
        rr1  = surface_list%flux_surfaces(i)%s(1,j);   ss1  = surface_list%flux_surfaces(i)%t(1,j)
        drr1 = surface_list%flux_surfaces(i)%s(2,j);   dss1 = surface_list%flux_surfaces(i)%t(2,j)
        rr2  = surface_list%flux_surfaces(i)%s(3,j);   ss2  = surface_list%flux_surfaces(i)%t(3,j)
        drr2 = surface_list%flux_surfaces(i)%s(4,j);   dss2 = surface_list%flux_surfaces(i)%t(4,j)
        i_elm = surface_list%flux_surfaces(i)%elm(j)

        call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        call interp_RZ(node_list,element_list,i_elm,rr2,ss2,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                            ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)
        dRRg1_dt = dRRg1_dr * drr1 + dRRg1_ds * dss1
        dZZg1_dt = dZZg1_dr * drr1 + dZZg1_ds * dss1
        dRRg2_dt = dRRg2_dr * drr2 + dRRg2_ds * dss2
        dZZg2_dt = dZZg2_dr * drr2 + dZZg2_ds * dss2
          
        do k=1,n_sub
          rr = 2.d0 * real(k-1)/real(n_sub-1) - 1.d0
          call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, rr, R, dRR_dr)
          call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, rr, Z, dZZ_dr)
          
          write(101,'(A,i2,A,f15.4)')         ' rplot[',k-1,'] = ',R
          write(101,'(A,i2,A,f15.4)')         ' zplot[',k-1,'] = ',Z
        enddo
        write(101,'(A)')                      ' pylab.plot(rplot,zplot, "r")'
      enddo
    enddo
  close(101)



end subroutine print_py_plot_unordered_flux_surfaces















! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface using the fact that they are ordered (ie. from the first piece of a part until its last one)
subroutine print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list, colour, dashed)

  use data_structure
  use mod_interp, only: interp_RZ
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)          :: filename
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  type (type_surface_list), intent(in)          :: surface_list
  character*1,              intent(in)          :: colour
  logical,                  intent(in)          :: dashed
  
  ! --- Internal variables
  integer       :: i, j, k, l, n_sub, count
  integer       :: i_elm
  real*8        :: st, rr, ss, dr_flux, ds_flux
  real*8        :: R, Z
  real*8        :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2
  
  n_sub = 2 ! minimum 2
  
  open(101,file=filename,position='append')
    write(101,'(A,i6,A)')' r = N.zeros(',n_pieces_max,')'
    write(101,'(A,i6,A)')' z = N.zeros(',n_pieces_max,')'
    do i=1,surface_list%n_psi
      do j=1,surface_list%flux_surfaces(i)%n_parts
        count = 0
        do k = surface_list%flux_surfaces(i)%parts_index(j), surface_list%flux_surfaces(i)%parts_index(j+1)-1
          do l = 1,n_sub
            rr1  = surface_list%flux_surfaces(i)%s(1,k);   ss1  = surface_list%flux_surfaces(i)%t(1,k)
            drr1 = surface_list%flux_surfaces(i)%s(2,k);   dss1 = surface_list%flux_surfaces(i)%t(2,k)
            rr2  = surface_list%flux_surfaces(i)%s(3,k);   ss2  = surface_list%flux_surfaces(i)%t(3,k)
            drr2 = surface_list%flux_surfaces(i)%s(4,k);   dss2 = surface_list%flux_surfaces(i)%t(4,k)
            st = -1.d0 + 2.d0*real(l-1)/real(n_sub-1)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            i_elm = surface_list%flux_surfaces(i)%elm(k)
            call interp_RZ(node_list,element_list,i_elm,rr,ss,R,Z)
            write(101,'(A,i6,A,f15.4)')' r[',count,'] = ',R
            write(101,'(A,i6,A,f15.4)')' z[',count,'] = ',Z
            count = count + 1
          enddo
        enddo
        write(101,'(A,i6)')' n_points = ', count
        if (dashed) then
          write(101,'(A,A,A)')' pylab.plot(r[0:n_points],z[0:n_points], "',colour,'--")'
        else
          write(101,'(A,A,A)')' pylab.plot(r[0:n_points],z[0:n_points], "',colour,'")'
        endif
        write(101,'(A,A,A)')' pylab.plot(r[0],z[0], "',colour,'x")'
        write(101,'(A,A,A)')' pylab.plot(r[n_points-1],z[n_points-1], "',colour,'x")'
      enddo
    enddo
  close(101)

end subroutine print_py_plot_ordered_flux_surfaces









  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots the wall
subroutine print_py_plot_wall(filename)

  use data_structure
  use phys_module
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)          :: filename
  
  ! --- Internal variables
  integer       :: i
  
  open(101,file=filename,position='append')
    write  (101,'(A,i8,A)')                                             ' r_points = N.zeros(',n_limiter,')'
    write  (101,'(A,i8,A)')                                             ' z_points = N.zeros(',n_limiter,')'
    do i=1,n_limiter
      write(101,'(A,i8,A,f15.4)')                                       ' r_points[',i-1,'] = ',R_limiter(i)
      write(101,'(A,i8,A,f15.4)')                                       ' z_points[',i-1,'] = ',Z_limiter(i)
    enddo
    write(101,'(A,A,A)')                                                ' pylab.plot(r_points,z_points, "b")'
  close(101)

end subroutine print_py_plot_wall













  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface using the fact that they are ordered (ie. from the first piece of a part until its last one)
subroutine print_py_plot_points(filename, n_points, R_points, Z_points)

  use data_structure
  use phys_module
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)          :: filename
  integer,                  intent(in)          :: n_points
  real*8,                   intent(in)          :: R_points(n_points), Z_points(n_points)
  
  ! --- Internal variables
  integer       :: i
  
  open(101,file=filename,position='append')
    do i=1,n_points
      write(101,'(A,f15.4)')                                            ' r_points = ',R_points(i)
      write(101,'(A,f15.4)')                                            ' z_points = ',Z_points(i)
      write(101,'(A)')                                                  ' pylab.plot(r_points,z_points, "xk", markersize=10)'
    enddo
  close(101)

end subroutine print_py_plot_points















  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface using the fact that they are ordered (ie. from the first piece of a part until its last one)
subroutine print_py_plot_line(filename, n_points, R_points, Z_points, colour, dashed)

  use data_structure
  use phys_module
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)          :: filename
  integer,                  intent(in)          :: n_points
  real*8,                   intent(in)          :: R_points(n_points), Z_points(n_points)
  character*1,              intent(in)          :: colour
  logical,                  intent(in)          :: dashed
  
  ! --- Internal variables
  integer       :: i
  
  open(101,file=filename,position='append')
    write  (101,'(A,i8,A)')                                             ' r_points = N.zeros(',n_points,')'
    write  (101,'(A,i8,A)')                                             ' z_points = N.zeros(',n_points,')'
    do i=1,n_points
      write(101,'(A,i8,A,f17.9)')                                       ' r_points[',i-1,'] = ',R_points(i)
      write(101,'(A,i8,A,f17.9)')                                       ' z_points[',i-1,'] = ',Z_points(i)
    enddo
    if (dashed) then
      write(101,'(A,A,A)')                                              ' pylab.plot(r_points,z_points, "',colour,'--")'
    else
      write(101,'(A,A,A)')                                              ' pylab.plot(r_points,z_points, "',colour,'")'
    endif
  close(101)

end subroutine print_py_plot_line


end module py_plots_grids

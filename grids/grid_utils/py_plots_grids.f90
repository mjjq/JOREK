module py_plots_grids
  implicit none
contains

!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine print_py_plot_prepare_plot(filename)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)		:: filename
  
  open(101,file=filename)
    write(101,'(A)')			      '#!/usr/bin/env python'
    write(101,'(A)')			      'import numpy as N'
    write(101,'(A)')			      'import pylab'
    write(101,'(A)')			      'def main():'
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
  character*256,            intent(in)		:: filename
  
  open(101,file=filename,position='append')
    write(101,'(A)')			      ' pylab.axis("equal")'
    write(101,'(A)')			      ' pylab.show()'
    write(101,'(A)')			      ' '
    write(101,'(A)')			      'main()'
  close(101)



end subroutine print_py_plot_finish_plot






  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface not using the fact that they are ordered (ie. plotting piece after piece)
subroutine print_py_plot_unordered_flux_surfaces(filename, node_list, element_list, surface_list)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)		:: filename
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface_list), intent(inout)	:: surface_list
  
  ! --- Internal variables
  integer	:: i, j
  integer	:: i_elm
  real*8	:: rr,    ss
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  
  open(101,file=filename,position='append')
    write(101,'(A)')			      ' rplot = N.zeros(2)'
    write(101,'(A)')			      ' zplot = N.zeros(2)'
    do i=1,surface_list%n_psi
      do j=1,surface_list%flux_surfaces(i)%n_pieces
  	i_elm = surface_list%flux_surfaces(i)%elm(j)
  	rr    = surface_list%flux_surfaces(i)%s(1,j)
  	ss    = surface_list%flux_surfaces(i)%t(1,j)
  	call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
        						  Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
  	write(101,'(A,f15.4)')  	      ' rplot[0] = ',R
  	write(101,'(A,f15.4)')  	      ' zplot[0] = ',Z
        i_elm = surface_list%flux_surfaces(i)%elm(j)
  	rr    = surface_list%flux_surfaces(i)%s(3,j)
  	ss    = surface_list%flux_surfaces(i)%t(3,j)
  	call interp_RZ(node_list,element_list,i_elm,rr,ss,R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
        						  Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
  	write(101,'(A,f15.4)')  	      ' rplot[1] = ',R
  	write(101,'(A,f15.4)')  	      ' zplot[1] = ',Z
  	write(101,'(A)')		      ' pylab.plot(rplot,zplot, "r")'
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
subroutine print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)		:: filename
  type (type_node_list),    intent(in)		:: node_list
  type (type_element_list), intent(in)		:: element_list
  type (type_surface_list), intent(in)		:: surface_list
  
  ! --- Internal variables
  integer	:: i, j, k
  integer	:: i_elm
  real*8	:: rr,    ss
  real*8	:: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8	:: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  
  open(101,file=filename,position='append')
    write(101,'(A,i6,A)')						' r = N.zeros(',n_pieces_max,')'
    write(101,'(A,i6,A)')						' z = N.zeros(',n_pieces_max,')'
    do i=1,surface_list%n_psi
      do j=1,surface_list%flux_surfaces(i)%n_parts
        write(101,'(A,i6)')						' n_points = ', surface_list%flux_surfaces(i)%parts_index(j+1) &
	                                                                               -surface_list%flux_surfaces(i)%parts_index(j  ) + 1
        do k = surface_list%flux_surfaces(i)%parts_index(j), surface_list%flux_surfaces(i)%parts_index(j+1)-1
    	  rr	= surface_list%flux_surfaces(i)%s(1,k)
    	  ss	= surface_list%flux_surfaces(i)%t(1,k)
    	  i_elm = surface_list%flux_surfaces(i)%elm(k)
    	  call interp_RZ(node_list,element_list,i_elm,rr,ss,&
	                 R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
  			 Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
    	  write(101,'(A,i6,A,f15.4)')					' r[',k-surface_list%flux_surfaces(i)%parts_index(j),'] = ',R
    	  write(101,'(A,i6,A,f15.4)')					' z[',k-surface_list%flux_surfaces(i)%parts_index(j),'] = ',Z
        enddo
    	rr    = surface_list%flux_surfaces(i)%s(3,surface_list%flux_surfaces(i)%parts_index(j+1)-1)
    	ss    = surface_list%flux_surfaces(i)%t(3,surface_list%flux_surfaces(i)%parts_index(j+1)-1)
    	i_elm = surface_list%flux_surfaces(i)%elm(surface_list%flux_surfaces(i)%parts_index(j+1)-1)
    	call interp_RZ(node_list,element_list,i_elm,rr,ss,&
	  	       R,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
  	               Z,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
    	write(101,'(A,i6,A,f15.4)')					' r[',surface_list%flux_surfaces(i)%parts_index(j+1)-surface_list%flux_surfaces(i)%parts_index(j),'] = ',R
    	write(101,'(A,i6,A,f15.4)')					' z[',surface_list%flux_surfaces(i)%parts_index(j+1)-surface_list%flux_surfaces(i)%parts_index(j),'] = ',Z
        write(101,'(A)')						' pylab.plot(r[0:n_points],z[0:n_points], "r")'
      enddo
    enddo
  close(101)

end subroutine print_py_plot_ordered_flux_surfaces









  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This plots fluxsurface using the fact that they are ordered (ie. from the first piece of a part until its last one)
subroutine print_py_plot_wall(filename)

  use data_structure
  use phys_module
  implicit none
  
  ! --- Routine parameters
  character*256,            intent(in)		:: filename
  
  ! --- Internal variables
  integer	:: i
  
  open(101,file=filename,position='append')
    write(101,'(A)')							' r_wall = N.zeros(2)'
    write(101,'(A)')							' z_wall = N.zeros(2)'
    if (n_limiter .ne. 0) then
      do i=1,n_limiter-1
    	write(101,'(A,f15.4)')						' r_wall[0] = ',R_limiter(i)
    	write(101,'(A,f15.4)')						' z_wall[0] = ',Z_limiter(i)
    	write(101,'(A,f15.4)')						' r_wall[1] = ',R_limiter(i+1)
    	write(101,'(A,f15.4)')						' z_wall[1] = ',Z_limiter(i+1)
        write(101,'(A)')						' pylab.plot(r_wall,z_wall, "b")'
      enddo
    endif
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
  character*256,            intent(in)		:: filename
  integer,                  intent(in)		:: n_points
  real*8,                   intent(in)		:: R_points(n_points), Z_points(n_points)
  
  ! --- Internal variables
  integer	:: i
  
  open(101,file=filename,position='append')
    do i=1,n_points
      write(101,'(A,f15.4)')						' r_points = ',R_points(i)
      write(101,'(A,f15.4)')						' z_points = ',Z_points(i)
      write(101,'(A)')							' pylab.plot(r_points,z_points, "xk", markersize=10)'
    enddo
  close(101)

end subroutine print_py_plot_points


end module py_plots_grids

!> The routine finds fluxsurfaces by finding crossings with the edges of the elements
subroutine find_flux_surfaces(xpoint,xcase,node_list,element_list,surface_list)

use constants
use tr_module 
use data_structure

implicit none

! --- Routine parameters
logical,                  intent(in)     :: xpoint
integer,                  intent(in)     :: xcase
type (type_node_list)   , intent(in)	 :: node_list
type (type_element_list), intent(in)	 :: element_list
type (type_surface_list), intent(inout)  :: surface_list

! --- Local variables
real*8  :: psimin, psimax, a0, a1, a2, a3
real*8  :: dpsi_dr(4),dpsi_ds(4)
real*8  :: p1, dp1, dp4, p4, p2, p3, r_psi(4), s_psi(4), tht(4)
real*8  :: s, s2, s3, r_tmp, s_tmp, psr_tmp, pss_tmp, ttmp, tt
real*8  :: psi_xpoint(2),R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2), r_av, s_av

real*8  :: RRg(4),dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg(4),dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
integer :: l, i_neigh, Xneigh, icount
integer :: my_id, i, j, k, ifound, iv, im, is, n1, n2, n3
integer :: ifail, itht(4), itmp,i_elm_xpoint(2)

write(*,*) '***********************************'
write(*,*) '*   find_flux_surfaces            *'
write(*,*) '***********************************'
!write(*,*) ' n_psi : ',surface_list%n_psi
!write(*,*) ' values : ',surface_list%psi_values(1),surface_list%psi_values(surface_list%n_psi)

my_id = 1 ! Just don't want the printout...

if (allocated(surface_list%flux_surfaces)) then
   call tr_unregister_mem(sizeof(surface_list%flux_surfaces),"surface_list%flux_surfaces")
   deallocate(surface_list%flux_surfaces)
end if

allocate(surface_list%flux_surfaces(1:surface_list%n_psi))
call tr_register_mem(sizeof(surface_list%flux_surfaces),"surface_list%flux_surfaces")

do j=1, surface_list%n_psi
  surface_list%flux_surfaces(j)%n_pieces = 0
  surface_list%flux_surfaces(j)%elm      = 0
  surface_list%flux_surfaces(j)%s        = 0
  surface_list%flux_surfaces(j)%t        = 0
enddo

if (xpoint) call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)



do i=1, element_list%n_elements
       
  call psi_minmax(node_list,element_list,i,psimin,psimax)
 
  do j=1, surface_list%n_psi

    ifound = 0

!    if ((surface_list%psi_values(j) .ge. psimin) .and. (surface_list%psi_values(j) .le. psimax)) then

      do iv=1, 4

        im = MOD(iv,4) + 1
        n1 = element_list%element(i)%vertex(iv)
        n2 = element_list%element(i)%vertex(im)

        is = mod(iv+1,2) + 2

        p1  =  node_list%node(n1)%values(1,1,1)  * element_list%element(i)%size(iv,1)
        dp1 =  node_list%node(n1)%values(1,is,1) * element_list%element(i)%size(iv,is)
        p4  =  node_list%node(n2)%values(1,1,1)  * element_list%element(i)%size(im,1)
        dp4 =  node_list%node(n2)%values(1,is,1) * element_list%element(i)%size(im,is)

        p2  = p1 + dp1
        p3  = p4 + dp4

        a3 = -        p1 + 3.d0 * p2 - 3.d0 * p3 + p4
        a2 = + 3.d0 * p1 - 6.d0 * p2 + 3.d0 * p3
        a1 = - 3.d0 * p1 + 3.d0 * p2
        a0 =          p1                                       - surface_list%psi_values(j)

        call SOLVP3(a0,a1,a2,a3,s,s2,s3,ifail)
          
        if ((s .ge. 0.d0) .and. (s .le. 1.d0)) then

          ifound = ifound + 1

!         write(*,*) ' first solution : ',s

          call flux_surface_add_point(node_list,element_list,surface_list,s,i,iv,ifound,r_psi,s_psi,dpsi_dr,dpsi_ds)

        endif

        if ((s2 .ge. 0.d0) .and. (s2 .le. 1.d0)) then

          ifound = ifound + 1

!         write(*,*) ' second solution : ',s2

          call flux_surface_add_point(node_list,element_list,surface_list,s2,i,iv,ifound,r_psi,s_psi,dpsi_dr,dpsi_ds)

          if (abs(s3) .le. 1.d0)         write(*,*) ' WARNING another solution : ',s3

        endif

      enddo ! end of 4 edges
  
      if (ifound .eq. 2) then

        call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(1:2),s_psi(1:2),dpsi_dr(1:2),dpsi_ds(1:2))
	
      elseif (ifound .eq. 4) then
      
! complicated : 2 line pieces but which point belongs to which line piece?

        r_av = (r_psi(1)+r_psi(2)+r_psi(3)+r_psi(4))/4.d0
        s_av = (s_psi(1)+s_psi(2)+s_psi(3)+s_psi(4))/4.d0

        tht(1:4) = atan2(s_psi(1:4)-s_av,r_psi(1:4)-r_av)

        where (tht .lt. 0.d0) tht = tht + 2.d0*PI

        itht(1)= 1; itht(2) = 2; itht(3) = 3; itht(4) = 4
	
	if (tht(2) .lt. tht(1)) then
	  itmp = itht(1); itht(1) = itht(2) ; itht(2) = itmp;
	endif
	if (tht(4) .lt. tht(3)) then
	  itmp = itht(3); itht(3) = itht(4) ; itht(4) = itmp;
	endif
	if (tht(itht(3)) .lt. tht(itht(2))) then
	  itmp = itht(2); itht(2) = itht(3) ; itht(3) = itmp;
	endif
	if (tht(itht(2)) .lt. tht(itht(1))) then
	  itmp = itht(1); itht(1) = itht(2) ; itht(2) = itmp;
	endif
	if (tht(itht(4)) .lt. tht(itht(3))) then
	  itmp = itht(3); itht(3) = itht(4) ; itht(4) = itmp;
	endif
        if (tht(itht(3)) .lt. tht(itht(2))) then
          itmp = itht(2); itht(2) = itht(3) ; itht(3) = itmp;
        endif
        
	
        if ((xpoint) .and. (     ((i .eq. i_elm_xpoint(1)) .and. (xcase .ne. 2) .and. (surface_list%psi_values(j) .eq. psi_xpoint(1)) )  &
	                    .or. ((i .eq. i_elm_xpoint(2)) .and. (xcase .ne. 1) .and. (surface_list%psi_values(j) .eq. psi_xpoint(2)) )  ) ) then

	  call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(1:3:2)), &
          	   s_psi(itht(1:3:2)),dpsi_dr(itht(1:3:2)),dpsi_ds(itht(1:3:2)))
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(2:4:2)), &
                                s_psi(itht(2:4:2)),dpsi_dr(itht(2:4:2)),dpsi_ds(itht(2:4:2)))

        else

          ! This is a little tricky, we look if the element is a neighboor of one of the Xpoints
	  Xneigh = 0
	  do k=1,4
	    i_neigh = element_list%element(i)%neighbours(k)
	    if( (xcase .ne. 2) .and. (i_neigh .eq. i_elm_xpoint(1)) ) then
	      Xneigh = 1
	      exit
	    endif
	    if( (xcase .ne. 1) .and. (i_neigh .eq. i_elm_xpoint(2)) ) then
	      Xneigh = 2
	      exit
	    endif
	  enddo
          ! If it is a neighboor, then record all four intersections (also do that for cases where
	  ! the element is i_elm_xpoint, but the flux surface is not the LCFS)
	  if( (Xneigh .gt. 0) &
	    .or. ((i .eq. i_elm_xpoint(1)) .and. (xcase .ne. 2)) & 
	    .or. ((i .eq. i_elm_xpoint(2)) .and. (xcase .ne. 1)) ) then
	    do k=1,4
	      call interp_RZ(node_list,element_list,i,r_psi(k),s_psi(k),&
	     		     RRg(k),dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,	    &
	       		     ZZg(k),dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
	    enddo
	  endif
          ! Then, look if the element is above/below or right/left of i_elm_xpoint, 
	  ! and then reorder the points 1,2,3,4 so that 1,2 are always right/above Xpoint,
	  ! and 3,4 are always left/below Xpoint
	  if(Xneigh .gt. 0) then
	    if( (maxval(RRg) .gt. R_xpoint(Xneigh)) .and. (minval(RRg) .lt. R_xpoint(Xneigh)) ) then
	      icount = 0
	      do k=1,4
	        if(RRg(k) .gt. R_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	      do k=1,4
	        if(RRg(k) .lt. R_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	    else
	      icount = 0
	      do k=1,4
	        if(ZZg(k) .gt. Z_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	      do k=1,4
	        if(ZZg(k) .lt. Z_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	    endif
	  endif
          
	  ! In the case where the element actually is i_elm_xpoint, 
	  ! but the flux surface is not the LCFS, we need to check if the line is right&left
	  ! or above&below the Xpoint
	  if( (Xneigh .eq. 0) &
	    .and. (    ((i .eq. i_elm_xpoint(1)) .and. (xcase .ne. 2)) &
	          .or. ((i .eq. i_elm_xpoint(2)) .and. (xcase .ne. 1)) ) ) then
	    if(i .eq. i_elm_xpoint(1)) Xneigh = 1
	    if(i .eq. i_elm_xpoint(2)) Xneigh = 2
	    if(surface_list%psi_values(j) .gt. psi_xpoint(Xneigh)) then
	      icount = 0
	      do k=1,4
	        if(RRg(k) .gt. R_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	      do k=1,4
	        if(RRg(k) .lt. R_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	    else
	      icount = 0
	      do k=1,4
	        if(ZZg(k) .gt. Z_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	      do k=1,4
	        if(ZZg(k) .lt. Z_xpoint(Xneigh)) then
	          icount = icount + 1
	          itht(icount) = k
	        endif
	      enddo
	    endif
	  endif

          ! Then add the lines
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(1:2)), &
                              s_psi(itht(1:2)),dpsi_dr(itht(1:2)),dpsi_ds(itht(1:2)))
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(3:4)), &
                              s_psi(itht(3:4)),dpsi_dr(itht(3:4)),dpsi_ds(itht(3:4)))

        endif

      endif
     
!    endif

  enddo

enddo

return
end

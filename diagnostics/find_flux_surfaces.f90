subroutine find_flux_surfaces(xpoint,node_list,element_list,surface_list)
!-----------------------------------------------------------------------
! finds fluxsurfaces by finding crossings at the edge of an element
!-----------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

real*8  :: psimin, psimax, a0, a1, a2, a3, PI
real*8  :: psi_test, dpsi_dr(4),dpsi_ds(4), dRR_dr(4), dRR_ds(4), dZZ_dr(4), dZZ_ds(4)
real*8  :: p1, dp1, dp4, p4, p2, p3, RR_psi(4), ZZ_psi(4), r_psi(4), s_psi(4), tht(4)
real*8  :: s, s2, s3, r_tmp, s_tmp, psr_tmp, pss_tmp, ttmp, tt
real*8  :: psi_xpoint,R_xpoint,Z_xpoint,s_xpoint,t_xpoint, r_av, s_av
real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss, RZ_jac, PSI_R, pSI_Z

integer :: i, j, k, ifound, iv, im, is, n1, n2, n3
integer :: ifail, itht(4), itmp,i_elm_xpoint
logical :: xpoint

write(*,*) '***********************************'
write(*,*) '*   find_flux_surfaces            *'
write(*,*) '***********************************'
write(*,*) ' n_psi : ',surface_list%n_psi

PI = 2.d0* asin(1.d0)

if (allocated(surface_list%flux_surfaces)) deallocate(surface_list%flux_surfaces)

allocate(surface_list%flux_surfaces(surface_list%n_psi))
!print*,"findflux deb"
do j=1, surface_list%n_psi
  surface_list%flux_surfaces(j)%n_pieces = 0
  surface_list%flux_surfaces(j)%elm      = 0
  surface_list%flux_surfaces(j)%s        = 0
  surface_list%flux_surfaces(j)%t        = 0
!print*,"findflux fin",j,surface_list%n_psi
enddo

if (xpoint) call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)



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

!        write(*,*) ' found 4 points '

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

 !       write(*,*) i,j
 !       write(*,'(4f12.8)') tht
 !       write(*,'(4i5)') itht
 !       write(*,'(4f12.8)') tht(itht)
 !       write(*,'(4f12.8)') r_psi
 !       write(*,'(4f12.8)') s_psi
 !       write(*,'(4f12.8)') dpsi_dr
 !       write(*,'(4f12.8)') dpsi_ds

        if ((xpoint) .and. (i .eq. i_elm_xpoint)) then
          write(*,*) ' adding x-point'

          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(1:3:2)), &
                                s_psi(itht(1:3:2)),dpsi_dr(itht(1:3:2)),dpsi_ds(itht(1:3:2)))
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(2:4:2)), &
                                s_psi(itht(2:4:2)),dpsi_dr(itht(2:4:2)),dpsi_ds(itht(2:4:2)))

          do k=1,4
            call interp_RZ(node_list,element_list,i_elm_xpoint,r_psi(k),s_psi(k), &
                           RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                           ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

            call interp(node_list,element_list,i_elm_xpoint,1,1,r_psi(k),s_psi(k),&
                        PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

            RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

            PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
            PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

            write(*,*) 'x-point : ',r_psi(k),s_psi(k),PSI_R,PSI_Z

          enddo


        else

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

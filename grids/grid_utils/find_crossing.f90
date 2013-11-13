subroutine find_crossing(node_list,element_list,surface_list,j_surf,R_c,Z_c, &
                         R_out,Z_out,ielm_flux,r_flux,s_flux,t_tht,ifail)
!-------------------------------------------------------------------------
! solves two non-linear equations using Newtons method (from numerical recipes)
! LU decomposition replaced by explicit solution of 2x2 matrix.
!
! finds the crossing of two coordinate lines given as a series of cubics
!-------------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer :: i, k, j_surf, ifail, ntrial, istart
integer :: k_flux, i_elm, ielm_flux
real*8  :: R_out,Z_out,r_flux,s_flux,t_flux,t_tht,dr_flux, ds_flux
real*8  :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2
real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8  :: RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss
real*8  :: ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss
real*8  :: dRRg1_dt, dZZg1_dt, dRRg2_dt, dZZg2_dt
real*8  :: RR_flux, dRR_flux, RR_tht, dRR_tht,  ZZ_flux, dZZ_flux, ZZ_tht, dZZ_tht
real*8  :: R_c(4), Z_c(4), x(2), FVEC(2), FJAC(2,2), p(2), x_previous(2)
real*8  :: tolx, tolf, errx, errf, temp, dis, max_step

if ((R_c(1) .eq. R_c(3)) .and. (Z_c(1) .eq. Z_c(3))) then
  ifail = 9
  return
endif

ielm_flux = 0

do k=1,surface_list%flux_surfaces(j_surf)%n_pieces

  rr1  = surface_list%flux_surfaces(j_surf)%s(1,k);   ss1  = surface_list%flux_surfaces(j_surf)%t(1,k)
  drr1 = surface_list%flux_surfaces(j_surf)%s(2,k);   dss1 = surface_list%flux_surfaces(j_surf)%t(2,k)
  rr2  = surface_list%flux_surfaces(j_surf)%s(3,k);   ss2  = surface_list%flux_surfaces(j_surf)%t(3,k)
  drr2 = surface_list%flux_surfaces(j_surf)%s(4,k);   dss2 = surface_list%flux_surfaces(j_surf)%t(4,k)

  i_elm = surface_list%flux_surfaces(j_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
  call interp_RZ(node_list,element_list,i_elm,rr2,ss2,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                      ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)
  dRRg1_dt = dRRg1_dr * drr1 + dRRg1_ds * dss1
  dZZg1_dt = dZZg1_dr * drr1 + dZZg1_ds * dss1
  dRRg2_dt = dRRg2_dr * drr2 + dRRg2_ds * dss2
  dZZg2_dt = dZZg2_dr * drr2 + dZZg2_ds * dss2

  ntrial = 20
  tolx = 1.d-6
  tolf = 1.d-12

  do istart = 1,5

    if (istart .eq. 1) then
      x(1) = 0.d0
      x(2) = 0.d0
    elseif (istart .eq. 2) then
      x(1) = -0.71d0
      x(2) = -0.71d0
    elseif (istart .eq. 3) then
      x(1) =  0.71d0
      x(2) = -0.71d0
    elseif (istart .eq. 4) then
      x(1) =  0.71d0
      x(2) =  0.71d0
    elseif (istart .eq. 5) then
      x(1) = -0.71d0
      x(2) =  0.71d0
    endif

    ifail = 999
    max_step = 1.005d0

    do i=1,ntrial

      t_flux = x(1)
      t_tht  = x(2)

      call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, t_flux, RR_flux, dRR_flux)
      call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, t_flux, ZZ_flux, dZZ_flux)

      call CUB1D(R_c(1), R_c(2), R_c(3), R_c(4), t_tht, RR_tht, dRR_tht)
      call CUB1D(Z_c(1), Z_c(2), Z_c(3), Z_c(4), t_tht, ZZ_tht, dZZ_tht)

      FVEC(1)   = RR_tht - RR_flux
      FVEC(2)   = ZZ_tht - ZZ_flux
      FJAC(1,1) = - dRR_flux
      FJAC(1,2) =   dRR_tht
      FJAC(2,1) = - dZZ_flux
      FJAC(2,2) =   dZZ_tht

      errf=abs(fvec(1))+abs(fvec(2))

!    if (i .eq. ntrial) write(*,'(A,i3,4e16.8)') ' newton   : ',i,errf,errx,x

      if (errf .le. tolf) then

        t_flux = x(1)
        t_tht  = x(2)

        call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
        call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)

        k_flux    = k
        ielm_flux = i_elm
        R_out     = 0.5d0*(RR_tht + RR_flux)
        Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)

!     write(*,'(A,i3,4e16.8)') ' newton (1) : ',i,errf,errx,x

        ifail = 0
        return
      endif

      p = -fvec

      temp = p(1)
      dis  = fjac(2,2)*fjac(1,1)-fjac(1,2)*fjac(2,1)
      p(1) = (fjac(2,2)*p(1)-fjac(1,2)*p(2))/dis
      p(2) = (fjac(1,1)*p(2)-fjac(2,1)*temp)/dis

      errx=abs(p(1)) + abs(p(2))

      p = min(p,+0.25d0)
      p = max(p,-0.25d0)

      x = x + p

      ! Sometimes you need to look outside the grid...
      x = max(x,-max_step)
      x = min(x,+max_step)
      if(abs(x_previous(1)) .eq. max_step) max_step = max_step + 0.002d0
      if(abs(x_previous(2)) .eq. max_step) max_step = max_step + 0.002d0
      x_previous(1) = x(1)
      x_previous(2) = x(2)


      if (errx .le. tolx) then

        t_flux = x(1)
        t_tht  = x(2)

        call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
        call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)

        k_flux    = k
        ielm_flux = i_elm
        R_out     = 0.5d0*(RR_tht + RR_flux)
        Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)

!     write(*,'(A,i3,4e16.8)') ' newton (2) : ',i,errf,errx,x

        ifail = 0
        return
      endif

    enddo
  enddo

enddo

if (ielm_flux .eq. 0) ifail = 99

!write(*,'(A,8e16.8)') ' crossing wrong exit ',x,errx,errf

return
end

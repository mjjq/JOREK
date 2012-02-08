subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use constants
use data_structure
use phys_module
use mod_poiss
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm_axis, i_elm_xpoint(2), i_elm, ifail, xcase2
real*8     :: amplitude, psi, psi_axis, theta
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R_axis, Z_axis, s_axis, t_axis, R, Z, BigR, Ti0, Ti0_s, Te0, Te0_s, BigR_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi
real*8     :: zp, dp_dpsi, dp_dpsi2, dp_dz, dp_dz2, P_ss, P_st, P_tt, R_out,Z_out,s_out,t_out
real*8     :: psi_n, psi_bnd,psi_xpoint(2),R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2),psi_lim,R_lim,Z_lim
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t, ps0_x, ps0_y, R_s, R_t, Z_s, Z_t, xjac, direction, Btot
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions (400)       *'
  write(*,*) '***************************************'
endif

Z_xpoint(1) = -99.d0
Z_xpoint(2) = +99.d0

if (my_id .eq. 0) then

  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (ifail .ne. 0) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif

  psi_bnd = 0.d0
    
  if (xpoint2) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    psi_bnd  = psi_xpoint(1)
    if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
      psi_bnd = psi_xpoint(2)
    endif
    if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
    if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
  else
    Z_xpoint(1) = -999.d0
    Z_xpoint(2) = +999.d0
  endif

  if (freeboundary) then
    call find_limiter(node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
    if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
      psi_bnd = min(psi_lim,psi_bnd)
    endif
  endif

  if(xcase2 .ne. 2) write(*,'(A,3f10.5,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(1),ifail
  if(xcase2 .eq. 2) write(*,'(A,3f10.5,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(2),ifail

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1)
    Z   = node_list%node(i)%x(1,2)

    call density(      xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		       zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
    call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		       zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)

    call FFprime(      xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)

    zp       = zn * (zTi + zTe)
    dp_dpsi  = zn * (dTi_dpsi + dTe_dpsi) + dn_dpsi * (zTi + zTe)
    dp_dpsi2 = zn * (dTi_dpsi2 + dTe_dpsi2) + 2.d0 * dn_dpsi * (dTi_dpsi + dTe_dpsi) + dn_dpsi2 * (zTi + zTe)
    dp_dz    = zn * (dTi_dz + dTe_dz) + dn_dz * (zTi + zTe)
    dp_dz2   = zn * (dTi_dz2 + dTe_dz2) + 2.d0 * dn_dz * (dTi_dz + dTe_dz) + dn_dz2 * (zTi + zTe)				       

    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi  * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,5) = dn_dpsi  * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,5) = dn_dpsi  * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(4,2) &
                                    + dn_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

    node_list%node(i)%values(1,1,6) = zTi
    node_list%node(i)%values(1,2,6) = dTi_dpsi  * node_list%node(i)%values(1,2,1) + dTi_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,6) = dTi_dpsi  * node_list%node(i)%values(1,3,1) + dTi_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,6) = dTi_dpsi  * node_list%node(i)%values(1,4,1) + dTi_dz * node_list%node(i)%x(4,2) &
                                    + dTi_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dTi_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

    node_list%node(i)%values(1,1,8) = zTe
    node_list%node(i)%values(1,2,8) = dTe_dpsi  * node_list%node(i)%values(1,2,1) + dTe_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,8) = dTe_dpsi  * node_list%node(i)%values(1,3,1) + dTe_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,8) = dTe_dpsi  * node_list%node(i)%values(1,4,1) + dTe_dz * node_list%node(i)%x(4,2) &
                                    + dTe_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dTe_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

    node_list%node(i)%values(1,1,2) = - tauIC * zp 
    node_list%node(i)%values(1,2,2) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,2,1) + dp_dz * node_list%node(i)%x(2,2))
    node_list%node(i)%values(1,3,2) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,3,1) + dp_dz * node_list%node(i)%x(3,2))
    node_list%node(i)%values(1,4,2) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,4,1) + dp_dz * node_list%node(i)%x(4,2) &
                                    + dp_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dp_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2) )

    node_list%node(i)%values(1,1,7) = 0.d0        ! parallel velocity
    node_list%node(i)%values(1,2,7) = 0.d0
    node_list%node(i)%values(1,3,7) = 0.d0
    node_list%node(i)%values(1,4,7) = 0.d0
    
    node_list%node(i)%deltas = 0.d0
    
  enddo

endif


!----------------------------------------- flux boundary perturbation (to be completed, see Marina)
!if (my_id .eq. 0) then
!  do i=1,node_list%n_nodes
!    psi = node_list%node(i)%values(1,1,1)
!    R   = node_list%node(i)%x(1,1)
!    Z   = node_list%node(i)%x(1,2)
!    theta = atan2(Z-Z_axis, R-R_axis)
!    psi_bnd = 0.d0
!    if (xpoint2 .and. (xcase2 .ne. 2)) psi_bnd = psi_xpoint(1)
!    if (xpoint2 .and. (xcase2 .eq. 2)) psi_bnd = psi_xpoint(2)
!    if (xpoint2 .and. (xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) psi_bnd = psi_xpoint(2)
!    psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)     
!    if (node_list%node(i)%boundary .ne. 0) then
!      node_list%node(i)%values(2,1,1) =  0.01 * sin(2.d0*theta)
!      node_list%node(i)%values(3,1,1) = -0.01 * cos(2.d0*theta) 
!      node_list%node(i)%values(2,3,1) =  0.01 * cos(2.d0*theta) * 2.d0*PI/float(n_tht)
!      node_list%node(i)%values(3,3,1) = +0.01 * sin(2.d0*theta) * 2.d0*PI/float(n_tht)
!    endif
!  enddo
!call poisson(my_id,-2,node_list,element_list,1,3,2,xpoint2, xcase2) 
!call poisson(my_id,-2,node_list,element_list,1,3,3,xpoint2, xcase2) 
!endif    

!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,1)
      Z   = node_list%node(i)%x(1,2)
      psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)

      node_list%node(i)%values(in,1,4) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,4) = amplitude * (1. - 2.d0 * psi_n)/(psi_bnd - psi_axis) * node_list%node(i)%values(1,2,1)
      node_list%node(i)%values(in,3,4) = amplitude * (1. - 2.d0 * psi_n)/(psi_bnd - psi_axis) * node_list%node(i)%values(1,3,1)
      node_list%node(i)%values(in,4,4) = amplitude * (1. - 2.d0 * psi_n)/(psi_bnd - psi_axis) * node_list%node(i)%values(1,4,1)

      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. Z_xpoint(1)) .and. (xcase2 .ne. 2)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. Z_xpoint(2)) .and. (xcase2 .ne. 1)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               4,2,in, psi_axis,psi_bnd,xpoint2, xcase2,Z_xpoint,freeboundary,refinement,1)
enddo

!----------------------------------- fill in parallel velocity at boundary (on open field lines)
do i=1,node_list%n_nodes

  if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then

    ps0_s     = node_list%node(i)%values(1,2,1)
    ps0_t     = node_list%node(i)%values(1,3,1)
    R_s       = node_list%node(i)%x(2,1)
    R_t       = node_list%node(i)%x(3,1)
    Z_s       = node_list%node(i)%x(2,2)
    Z_t       = node_list%node(i)%x(3,2)

    xjac  =  R_s*Z_t - R_t*Z_s
    ps0_x = (	Z_t * ps0_s - Z_s * ps0_t ) / xjac
    ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

    direction = + ps0_x / abs(ps0_x)		 ! temporary solution for lower x-point only
    if (xcase2 .eq. 2) direction = -direction
    if ( (xcase2 .eq. 3) .and. (node_list%node(i)%x(1,2) .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) ) direction = -direction

    do in=1,n_tor

      BigR = node_list%node(i)%x(1,1)
      Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
      Ti0   = node_list%node(i)%values(in,1,6)
      Te0   = node_list%node(i)%values(in,1,8)
      node_list%node(i)%values(in,1,7) = direction / Btot * sqrt(GAMMA * (Ti0 + Te0))

      BigR_s = node_list%node(i)%x(2,1)
      Ti0_s   = node_list%node(i)%values(in,2,6)
      Te0_s   = node_list%node(i)%values(in,2,8)
      node_list%node(i)%values(in,2,7) = BigR_s / (BigR*Btot) * sqrt(GAMMA * (Ti0 + Te0)) + 0.5d0 / Btot * sqrt(GAMMA / (Ti0 + Te0)) * (Ti0_s + Te0_s)
      node_list%node(i)%values(in,2,7) = direction *  node_list%node(i)%values(in,2,7)

      if(xcase2 .eq. 1) then
        write(*,'(A,8e14.6)') ' Boundary condition (eq): ',BigR,psi_xpoint(1),node_list%node(i)%values(1,1,1),ps0_x,ps0_y, &
			    node_list%node(i)%values(in,1,n_var),BigR/F0 * sqrt(GAMMA*(Ti0 + Te0))
      endif
      if( (xcase2 .eq. 2) .or. ( (xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)) ) ) then
        write(*,'(A,8e14.6)') ' Boundary condition (eq): ',BigR,psi_xpoint(2),node_list%node(i)%values(1,1,1),ps0_x,ps0_y, &
			    node_list%node(i)%values(in,1,n_var),BigR/F0 * sqrt(GAMMA*(Ti0 + Te0))
      endif

    enddo
  endif
enddo

!call add_pellet(node_list,element_list,5.,0.08,0.03,R_geo-0.78,Z_geo)
!call add_pellet(node_list,element_list,25.,0.08,0.03,R_geo+0.85,Z_geo)

return
end

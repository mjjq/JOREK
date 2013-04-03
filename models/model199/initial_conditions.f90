subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use mod_poiss
implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_surface_list)     :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm_axis, i_elm_xpoint(2), ifail, i_elm, xcase2
real*8     :: amplitude, psi, psi_axis
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R_axis, Z_axis, s_axis, t_axis, R, Z, BigR
real*8     :: R_out, Z_out, s_out, t_out, R_lim, Z_lim, s_lim, t_lim, psi_lim
real*8     :: psi_n, psi_bnd,psi_xpoint(2),R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2)
real*8     :: p_s, p_t, p_ss, p_st, p_tt
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (199)      *'
  write(*,*) '***************************************'
endif

if (my_id .eq. 0) then

  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (ifail .ne. 0) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif

  psi_bnd  =   0.d0
  Z_xpoint(1) = -99.d0
  Z_xpoint(2) = +99.d0
    
  if (xpoint2) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    psi_bnd  = psi_xpoint(1)
    if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
      psi_bnd = psi_xpoint(2)
    endif
    if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
    if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
  endif

  if (freeboundary) then
    call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
    if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
      psi_bnd = min(psi_lim,psi_bnd)
    endif
  endif

  if(xcase2 .eq. 1) write(*,'(A,3f10.5,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(1),ifail
  if(xcase2 .eq. 2) write(*,'(A,3f10.5,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(2),ifail

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1)
    Z   = node_list%node(i)%x(1,2)

    call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
                                                               dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

    call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)

    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi  * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,5) = dn_dpsi  * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,5) = dn_dpsi  * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(4,2) &
                                    + dn_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

    node_list%node(i)%values(1,1,6) = zT
    node_list%node(i)%values(1,2,6) = dT_dpsi  * node_list%node(i)%values(1,2,1) + dT_dz * node_list%node(i)%x(2,2)
    node_list%node(i)%values(1,3,6) = dT_dpsi  * node_list%node(i)%values(1,3,1) + dT_dz * node_list%node(i)%x(3,2)
    node_list%node(i)%values(1,4,6) = dT_dpsi  * node_list%node(i)%values(1,4,1) + dT_dz * node_list%node(i)%x(4,2) &
                                    + dT_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dT_dz2   * node_list%node(i)%x(2,2)        * node_list%node(i)%x(3,2)

  enddo

endif

!---------------------------- initialise perturbations

amplitude = 1.d-10
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


!call add_pellet(node_list,element_list,50.,0.06,0.02,R_axis-0.96,Z_axis)

return
end

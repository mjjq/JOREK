subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2)
!-----------------------------------------------------------------------
! Solve the Grad-Shafranov equation to determine the plasma equilibrium
!   both freeboundary and fixed boundary solutions
!-----------------------------------------------------------------------
use tr_module 
use parameters
use data_structure
use phys_module
use mod_poiss
implicit none

          
! --- Routine parameters
integer,                      intent(in)    :: my_id
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_bnd_node_list),    intent(inout) :: bnd_node_list
type (type_bnd_element_list), intent(inout) :: bnd_elm_list
logical,                      intent(in)    :: xpoint2
integer,                      intent(in)    :: xcase2

! --- Local variables.
type (type_surface_list) :: surface_list, sep_list
integer    :: ierr, n_iter, iter, i, in, mm, i_elm_axis, i_elm_xpoint(2), i_elm_lim, ifail, i_elm
real*8     :: amplitude, psi, psi_bnd
real*8     :: zn,  dn_dpsi,  dn_dpsi2,  dn_dz,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi2_dz,  dn_dpsi_dz2
real*8     :: zT,  dT_dpsi,  dT_dpsi2,  dT_dz,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi2_dz,  dT_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: xx, x_s, x_t, x_st, x_ss, x_tt, yy, y_s, y_t, y_st, y_ss, y_tt
real*8     :: R_axis, Z_axis, s_axis, t_axis, psi_axis,R, Z, BigR, T0, BigR_s, T0_s
real*8     :: R_lim, Z_lim, s_lim, t_lim, psi_lim, R_out, Z_out, s_out, t_out
real*8     :: R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2), psi_xpoint(2)
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi, psi_n
real*8     :: ps0_s, ps0_t, p_s, p_t, p_ss, p_st, p_tt 
real*8     :: zj0_s, zj0_t, equil_error, equil_value, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
real*8     :: current_tot, current_ref, current_int, amix, ZKP, ZKI, ZKD, diff, R_xpoint2(2), Z_xpoint2(2)
real*8     :: sigmas(16)
integer    :: n_grids(10)
logical    :: freeboundary_equil2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*           equilibrium               *'
  write(*,*) '***************************************'
  write(*,*) '   freeboundary_equil : ',freeboundary_equil
  write(*,*) '   X-point      : ',xpoint2
  write(*,*) '   Xcase        : ',xcase2
endif

freeboundary_equil2 = freeboundary_equil
freeboundary_equil  = .false.

!------------------------------------ fixed boundary equilibrium
n_iter      = 200
psi_bnd     = 0.d0
Z_xpoint(1) = -99.d0
Z_xpoint(2) = +99.d0

do iter = 1, n_iter

  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if ((ifail .ne. 0) .and. (iter .le. 5)) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif
  
  if (xpoint2) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    if (ifail == 0) then ! (otherwise, keep the values of the previous iteration as a reasonable guess)
      psi_bnd  = psi_xpoint(1)
      if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
        psi_bnd = psi_xpoint(2)
      endif
      R_xpoint(1) = R_xpoint2(1)
      Z_xpoint(1) = Z_xpoint2(1)
      R_xpoint(2) = R_xpoint2(2)
      Z_xpoint(2) = Z_xpoint2(2)
      if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
      if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
    else
      if (freeboundary_equil) then
        Z_xpoint(1) = -99.d0
        Z_xpoint(2) = +99.d0
      endif
    endif
  endif

  call find_limiter(node_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
  if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
    if (psi_lim .lt. psi_bnd) then
      psi_bnd = psi_lim
      write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
    endif
  endif
    
  if(xcase2 .eq. 1) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(1),ifail
  if(xcase2 .eq. 2) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(2),ifail
  
  call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
               psi_axis,psi_bnd,xpoint2,xcase2,Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1

  diff = 0.d0
  do i=1, node_list%n_nodes
    diff = diff + abs(node_list%node(i)%deltas(1,1,1))
  enddo  
  diff = diff / float(node_list%n_nodes)
  
  write(*,'(A,I4,A,ES10.3)') ' Iteration ', iter, ': diff=', diff
  
  if ( (iter > 1) .and. (diff < 1.d-6) ) then
    write(*,'(A,I4,A)') ' Equilibrium converged: after', iter, ' iterations'
    exit
  else if ( iter == n_iter) then
    write(*,'(A,ES10.3)') ' Equilibrium not fully converged: diff=', diff
    exit
  end if

enddo


!--------------------------------------- freeboundary equilibrium
freeboundary_equil = freeboundary_equil2

current_int = 0.d0
amix        = 0.96
ZKP         = 1.0                ! PI feedback on the total current
ZKI         = 0.01 

if (freeboundary_equil) then

  n_iter =200
  
  call integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, xcase2, Z_xpoint, current_ref)
  
  do iter=1,n_iter
        
    call integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, xcase2, Z_xpoint, current_tot)
    
    current_int = current_int + (current_tot-current_ref)
    
    T_0  = T_0  * (1. - ZKP*(current_tot-current_ref)/current_ref - ZKI*current_int/current_ref)
    FF_0 = FF_0 * (1. - ZKP*(current_tot-current_ref)/current_ref - ZKI*current_int/current_ref)
    
    write(*,'(A,8e12.4)') 'FEEDBACK : ',current_ref,current_tot,current_int,T_0, FF_0
                   
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

    if ((ifail .ne. 0) .and. (iter .le. 5)) then
      call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
      call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
      write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
    endif
    
    psi_bnd = 0.d0
  
    if (xpoint2) then
      call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
      if (ifail .ne. 1) then      
        psi_bnd  = psi_xpoint(1)
        if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
          psi_bnd = psi_xpoint(2)
        endif
        if(xcase2 .eq. 1) Z_xpoint(2) = +99.d0
        if(xcase2 .eq. 2) Z_xpoint(1) = -99.d0
      else
        Z_xpoint(1) = -99.d0 
	Z_xpoint(2) = +99.d0
      endif
    endif

    call find_limiter(node_list,bnd_elm_list,psi_lim,R_lim,Z_lim)

    if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
      psi_bnd = min(psi_lim,psi_bnd)
      write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
    endif
       
    if(xcase2 .eq. 1) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(1),ifail
    if(xcase2 .eq. 2) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint(2),ifail

    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 psi_axis,psi_bnd,xpoint2,xcase2,Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1

!    call boundary_check
 
    diff = 0.d0
    do i=1, node_list%n_nodes
      diff = diff + abs(node_list%node(i)%deltas(1,1,1))
    enddo  
    diff = diff / float(node_list%n_nodes)
  
    write(*,'(A,i5,e14.6)') ' iteration, diff : ',iter,diff
  
    if ((iter .gt. 1) .and. (diff .lt. 1.d-6)) exit

  enddo

  if (my_id == 0) call boundary_check()

endif

!------------------------------- end of equilibrium, start filling data

do i=1,node_list%n_nodes

  psi = node_list%node(i)%values(1,1,1)
  R   = node_list%node(i)%x(1,1)
  Z   = node_list%node(i)%x(1,2)

  call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                             dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

  if (jorek_model .eq. 400) then
    call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
  		     zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)

    call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
  		     zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
    zT  	= zTi	       + zTe
    dT_dpsi	= dTi_dpsi     + dTe_dpsi
    dT_dpsi2	= dTi_dpsi2    + dTe_dpsi2
    dT_dpsi3	= dTi_dpsi3    + dTe_dpsi3
    dT_dz	= dTi_dz       + dTe_dz
    dT_dz2	= dTi_dz2      + dTe_dz2
    dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
    dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
    dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
  else
    call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
  		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
  endif

#ifdef fullmhd
  call F_profile(xpoint2, xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,&
                 F_prof,dF_dpsi,dF_dz, dF_dpsi2,dF_dz2,dF_dpsi_dz, &
                 zFFprime,dFFprime_dpsi,ddFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz)
		 
#else
  call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                             dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)
#endif

  zjz     = zFFprime      - R*R *      (dn_dpsi    * zT + zn * dT_dpsi)

  dj_dpsi = dFFprime_dpsi - R*R *      (dn_dpsi2   * zT + zn * dT_dpsi2  + 2.d0 * dn_dpsi * dT_dpsi)

  dj_dR   =               - 2.d0 * R * (dn_dpsi    * zT + zn * dT_dpsi)

  dj_dZ   = dFFprime_dz   - R*R *      (dn_dpsi_dz * zT + dn_dpsi * dT_dz + zn * dT_dpsi_dz + dn_dz * dT_dpsi)

  dj_dR_dR = - 2.d0     * (dn_dpsi     * zT + zn * dT_dpsi)

  dj_dZ_dZ = dFFprime_dz2   - R*R * ( dn_dpsi_dz2 * zT   + dn_dpsi_dz * dT_dz  + dn_dz * dT_dpsi_dz  + dn_dz2 * dT_dpsi &
                                    +  dn_dpsi_dz  * dT_dz + dn_dpsi    * dT_dz2 + zn    * dT_dpsi_dz2 + dn_dz  * dT_dpsi_dz)

  dj_dpsi2 = dFFprime_dpsi2 - R*R * (dn_dpsi3 * zT + 3.d0 * dn_dpsi * dT_dpsi2 + 3.d0 * dn_dpsi2 * dT_dpsi + zn * dT_dpsi3 )

  dj_dR_dZ   = - 2.d0 * R * (dn_dpsi_dz * zT + dn_dpsi * dT_dz + zn * dT_dpsi_dz + dn_dz * dT_dpsi)

  dj_dR_dpsi = - 2.d0 * R * (dn_dpsi2   * zT + zn * dT_dpsi2   + 2.d0 * dn_dpsi * dT_dpsi)

  dj_dZ_dpsi = dFFprime_dpsi_dz - R*R * ( dn_dpsi2_dz * zT    + dn_dz * dT_dpsi2     + 2.d0 * dn_dpsi_dz * dT_dpsi  &
                                          + dn_dpsi2    * dT_dz + zn    * dT_dpsi2_dz  + 2.d0 * dn_dpsi    * dT_dpsi_dz)

  node_list%node(i)%values(1,1,3) = zjz

  node_list%node(i)%values(1,2,3) = dj_dpsi * node_list%node(i)%values(1,2,1) &
                                  + dj_dR   * node_list%node(i)%x(2,1)        &
                                  + dj_dZ   * node_list%node(i)%x(2,2)

  node_list%node(i)%values(1,3,3) = dj_dpsi * node_list%node(i)%values(1,3,1) &
                                  + dj_dR   * node_list%node(i)%x(3,1)        &
                                  + dj_dZ   * node_list%node(i)%x(3,2)

  node_list%node(i)%values(1,4,3) = dj_dpsi  * node_list%node(i)%values(1,4,1) &
                                  + dj_dR    * node_list%node(i)%x(4,1)        &
                                  + dj_dZ    * node_list%node(i)%x(4,2)        &
                                  + dj_dR_dR * node_list%node(i)%x(2,1) * node_list%node(i)%x(3,1)  &
                                  + dj_dZ_dZ * node_list%node(i)%x(2,2) * node_list%node(i)%x(3,2)  &
                                  + dj_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                  + dj_dR_dZ * ( node_list%node(i)%x(2,1) * node_list%node(i)%x(3,2)          &
                                               + node_list%node(i)%x(3,1) * node_list%node(i)%x(2,2) )        &
                                  + dj_dR_dpsi*( node_list%node(i)%x(2,1) * node_list%node(i)%values(1,3,1)   &
                                               + node_list%node(i)%x(3,1) * node_list%node(i)%values(1,2,1) ) &
                                  + dj_dZ_dpsi*( node_list%node(i)%x(2,2) * node_list%node(i)%values(1,3,1)   &
                                               + node_list%node(i)%x(3,2) * node_list%node(i)%values(1,2,1) )

#ifdef fullmhd
  node_list%node(i)%psi_eq(:,:) = node_list%node(i)%values(1,:,:)
#endif

enddo

! --- Find flux surfaces and plot them; determine the q-profile.  
if (xpoint2 .and. (n_flux .gt. 1)) then
  
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
  
  n_grids = 0
  sigmas  = 0.d0
  
  ! Build up some arrays to send as routine parameters to define_flux_values
  sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
  sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
  sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
  sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
  sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
  sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
  sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv

  n_grids(1) = 2*n_flux   ; n_grids(2) = n_tht
  n_grids(3) = 2*n_open   ; n_grids(4) = 2*n_outer  ; n_grids(5) = 2*n_inner
  n_grids(6) = 2*n_private; n_grids(7) = 2*n_up_priv
  n_grids(8) = n_leg      ; n_grids(9) = n_up_leg
  if (xcase .eq. 1) then
    n_grids(4) = 0
    n_grids(5) = 0
    n_grids(7) = 0
    n_grids(9) = 0
  endif
  if (xcase .eq. 2) then
    n_grids(4) = 0
    n_grids(5) = 0
    n_grids(6) = 0
    n_grids(8) = 0
  endif

  ! Allocate surface_list structure (that's for plotting only)
  if (xcase2 .eq. 1) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_private
  if (xcase2 .eq. 2) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_up_priv
  if (xcase2 .eq. 3) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_outer + 2*n_inner + 2*n_private + 2*n_up_priv
  if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
  call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
  
  ! Allocate sep_list structure (that's for plotting only)  
  sep_list%n_psi =3
  if(xcase .eq. 3) sep_list%n_psi =6
  if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
  
  ! Define the flux values to be plotted...
  psi_axis = psi_axis+0.01 !Just offset a little, because finding surfaces along the side of an element (on the xpoint grid) can be hard...
  call define_flux_values(node_list, element_list, surface_list, sep_list, &
                          xcase2, R_xpoint, Z_xpoint, psi_xpoint, psi_axis, n_grids, sigmas)
  psi_axis = psi_axis-0.01 !Put it back, it's not used anyway, but just for principle!
  
else
  surface_list%n_psi = 200  
  if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
  call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
  
  do i = 1, surface_list%n_psi
    surface_list%psi_values(i) = 1.25d0*(float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis) + psi_axis
  enddo
  
  call find_flux_surfaces(xpoint2,xcase2,node_list,element_list,surface_list)

  sep_list%n_psi =1
  if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
  sep_list%psi_values(1) = psi_bnd

  call find_flux_surfaces(xpoint2,xcase2,node_list,element_list,sep_list)  
endif
  
if (freeboundary_equil) then
  call plot_coils(.true.)
  call plot_flux_surfaces(node_list,element_list,surface_list,.false.,4,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
  call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)

  call plot_flux_surfaces(node_list,element_list,surface_list,.true.,4,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
  call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
  call plot_coils(.false.)
else
  if (xpoint2 .and. (n_flux .gt. 1)) then
    call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
  else
    call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,.false.,0)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint2,xcase2)
  endif
endif
  
call q_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint)

if (allocated(surface_list%psi_values))    call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
if (allocated(surface_list%flux_surfaces)) deallocate(surface_list%flux_surfaces)
if (allocated(sep_list%psi_values))        call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)
 
return
end subroutine equilibrium

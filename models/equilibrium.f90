subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2)
!-----------------------------------------------------------------------
! Solve the Grad-Shafranov equation to determine the plasma equilibrium
!   both freeboundary and fixed boundary solutions
!-----------------------------------------------------------------------
use data_structure
use phys_module

implicit none

INTERFACE 
  SUBROUTINE POISSON(MY_ID,ITYPE,NODE_LIST,ELEMENT_LIST,           &
     &BND_NODE_LIST,BND_ELM_LIST,IVAR_IN,IVAR_OUT,I_HARM,PSI_AXIS, &
     &PSI_BND,XPOINT,Z_XPOINT,FREEBOUNDARY_EQUIL,REFINEMENT,ITER)
    USE DATA_STRUCTURE
    INTEGER(KIND=4), INTENT(IN) :: MY_ID
    INTEGER(KIND=4), INTENT(IN) :: ITYPE
    TYPE (TYPE_NODE_LIST), INTENT(INOUT) :: NODE_LIST
    TYPE (TYPE_ELEMENT_LIST), INTENT(INOUT) :: ELEMENT_LIST
    TYPE (TYPE_BND_NODE_LIST) :: BND_NODE_LIST
    TYPE (TYPE_BND_ELEMENT_LIST) :: BND_ELM_LIST
    INTEGER(KIND=4), INTENT(IN) :: IVAR_IN
    INTEGER(KIND=4), INTENT(IN) :: IVAR_OUT
    INTEGER(KIND=4), INTENT(IN) :: I_HARM
    REAL(KIND=8) :: PSI_AXIS
    REAL(KIND=8) :: PSI_BND
    LOGICAL(KIND=4), INTENT(IN) :: XPOINT
    REAL(KIND=8) :: Z_XPOINT
    LOGICAL(KIND=4), INTENT(IN) :: FREEBOUNDARY_EQUIL
    LOGICAL(KIND=4), INTENT(IN) :: REFINEMENT
    INTEGER(KIND=4), INTENT(IN) :: ITER
  END SUBROUTINE POISSON
END INTERFACE 
          
! --- Routine parameters
integer,                      intent(in)    :: my_id
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_bnd_node_list),    intent(inout) :: bnd_node_list
type (type_bnd_element_list), intent(inout) :: bnd_elm_list
logical,                      intent(in)    :: xpoint2

! --- Local variables.
type (type_surface_list) :: surface_list, sep_list
integer    :: ierr, n_iter, iter, i, in, mm, i_elm_axis, i_elm_xpoint, i_elm_lim, ifail, i_elm
real*8     :: amplitude, psi, psi_bnd
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: xx, x_s, x_t, x_st, x_ss, x_tt, yy, y_s, y_t, y_st, y_ss, y_tt
real*8     :: R_axis, Z_axis, s_axis, t_axis, psi_axis,R, Z, BigR, T0, BigR_s, T0_s
real*8     :: R_lim, Z_lim, s_lim, t_lim, psi_lim, R_out, Z_out, s_out, t_out
real*8     :: R_xpoint,Z_xpoint,s_xpoint,t_xpoint, psi_xpoint
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi, psi_n
real*8     :: ps0_s, ps0_t, p_s, p_t, p_ss, p_st, p_tt 
real*8     :: zj0_s, zj0_t, equil_error, equil_value, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
real*8     :: current_tot, current_ref, current_int, amix, ZKP, ZKI, ZKD, diff
logical    :: freeboundary_equil2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*           equilibrium               *'
  write(*,*) '***************************************'
  write(*,*) '   freeboundary_equil : ',freeboundary_equil
  write(*,*) '   X-point      : ',xpoint2
endif

freeboundary_equil2 = freeboundary_equil
freeboundary_equil  = .false.

!------------------------------------ fixed boundary equilibrium
n_iter   = 200
psi_bnd  = 0.d0
Z_xpoint = -99.d0

do iter = 1, n_iter

  call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if ((ifail .ne. 0) .and. (iter .le. 5)) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif
  
  if (xpoint2) then
    call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
    if (ifail == 0) then ! (otherwise, keep the values of the previous iteration as a reasonable guess)
      psi_bnd = psi_xpoint
    else
      if (freeboundary_equil) then
        Z_xpoint = -99.d0
      endif
    endif
  endif

  call find_limiter(node_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
  if (Z_lim .gt. Z_xpoint) then
    if (psi_lim .lt. psi_bnd) then
      psi_bnd = psi_lim
      write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
    endif
  endif
    
  write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint,ifail
  
  call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
               psi_axis,psi_bnd,xpoint2,Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1

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
  
  call integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, Z_xpoint, current_ref)
  
  do iter=1,n_iter
        
    call integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, Z_xpoint, current_tot)
    
    current_int = current_int + (current_tot-current_ref)
    
    T_0  = T_0  * (1. - ZKP*(current_tot-current_ref)/current_ref - ZKI*current_int/current_ref)
    FF_0 = FF_0 * (1. - ZKP*(current_tot-current_ref)/current_ref - ZKI*current_int/current_ref)
    
    write(*,'(A,8e12.4)') 'FEEDBACK : ',current_ref,current_tot,current_int,T_0, FF_0
                   
    call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

    if ((ifail .ne. 0) .and. (iter .le. 5)) then
      call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
      call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
      write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
    endif
    
    psi_bnd = 0.d0
  
    if (xpoint2) then
      call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
      if (ifail .ne. 1) then      
        psi_bnd = psi_xpoint
      else
        Z_xpoint = -99.d0
      endif
    endif

    call find_limiter(node_list,bnd_elm_list,psi_lim,R_lim,Z_lim)

    if (Z_lim .gt. Z_xpoint) then
      psi_bnd = min(psi_lim,psi_bnd)
      write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
    endif
       
    write(*,'(A,3e14.6,i3)') ' PSI_AXIS, PSI_BND : ',psi_axis,psi_bnd,Z_xpoint,ifail

    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 psi_axis,psi_bnd,xpoint2,Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1

!    call boundary_check
 
    diff = 0.d0
    do i=1, node_list%n_nodes
      diff = diff + node_list%node(i)%deltas(1,1,1)
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

  call density(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                             dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

  call temperature(xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
                                                             dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

  call FFprime(    xpoint2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                             dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz)
							      
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

enddo

! --- Find flux surfaces and plot them; determine the q-profile.

surface_list%n_psi = 200
  
if (allocated(surface_list%psi_values)) deallocate(surface_list%psi_values)
allocate(surface_list%psi_values(surface_list%n_psi))
  
if (xpoint2) then
  do i = 1, surface_list%n_psi
    surface_list%psi_values(i) =  1.5 * (float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis) + psi_axis
  enddo
else
  do i = 1, surface_list%n_psi
    surface_list%psi_values(i) = (float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis) + psi_axis
  enddo
endif
  
call find_flux_surfaces(xpoint2,node_list,element_list,surface_list)

sep_list%n_psi =1
if (allocated(sep_list%psi_values)) deallocate(sep_list%psi_values)
allocate(sep_list%psi_values(sep_list%n_psi))
sep_list%psi_values(1) = psi_bnd

call find_flux_surfaces(xpoint2,node_list,element_list,sep_list)  

if (freeboundary_equil) then
  call plot_coils(.true.)
  call plot_flux_surfaces(node_list,element_list,surface_list,.false.,4)
  call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1)

  call plot_flux_surfaces(node_list,element_list,surface_list,.true.,4)
  call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1)
  call plot_coils(.false.)
else
  call plot_flux_surfaces(node_list,element_list,surface_list,.true.,4)
  call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1)
endif
  
call q_profile(node_list,element_list,surface_list,psi_axis,psi_bnd,Z_xpoint)

if (allocated(surface_list%psi_values))    deallocate(surface_list%psi_values)
if (allocated(surface_list%flux_surfaces)) deallocate(surface_list%flux_surfaces)
if (allocated(sep_list%psi_values))        deallocate(sep_list%psi_values)
if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)
 
return
end subroutine equilibrium

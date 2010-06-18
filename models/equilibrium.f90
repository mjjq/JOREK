subroutine equilibrium(my_id,node_list,element_list,xpoint2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
implicit none

interface
   subroutine Poisson(my_id,itype,node_list,element_list,ivar_in,ivar_out,i_harm,xpoint)
     use data_structure

     type (type_node_list)    :: node_list
     type (type_element_list) :: element_list
     integer  :: my_id, itype, ivar_in, ivar_out, i_harm
     logical  :: xpoint
   end subroutine Poisson
end interface
type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer    :: ierr, my_id, n_iter, iter, i, in, mm, i_elm_axis, i_elm_xpoint
real*8     :: amplitude, psi, psi_axis
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: w, w_s, w_t, w_st, w_r, w_tht, w_r_tht, w_rr, w_tht_tht
real*8     :: xx, x_s, x_t, x_st, x_ss, x_tt, yy, y_s, y_t, y_st, y_ss, y_tt
real*8     :: rr, r2, r_s, r_t, r_st, tht, tht_s, tht_t, tht_st
real*8     :: r_x, r_y, r_xy, r_xx, r_yy, tht_x, tht_y, tht_xy, tht_xx, tht_yy
real*8     :: x_min, x_max, y_min, y_max,x ,y
real*8     :: R_axis, Z_axis, s_axis, t_axis, R, Z, BigR, T0, BigR_s, T0_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi, psi_n
real*8     :: psi_bnd,psi_xpoint,R_xpoint,Z_xpoint,s_xpoint,t_xpoint
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t, equil_error, equil_value, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*           equilibrium               *'
  write(*,*) '***************************************'
endif

n_iter = 51
i_elm_axis = 1 ! XL : uninitilised value... So let's say it'll be 1, to look in the first case of element array in interp.f90...

do iter=1,n_iter

  if (my_id .eq. 0) then

    call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis)
    if (xpoint2) call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)

  endif

  call poisson(my_id,-1,node_list,element_list,3,1,1,xpoint2)   !----------- for GS use -1

!  call poisson(my_id,2,node_list,element_list,3,1,1,xpoint2)   !----------- for cylinder use 2

enddo

if (my_id .eq. 0) then

  psi_bnd = 0.d0
  if (xpoint2) psi_bnd = psi_xpoint

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

endif

call poisson(my_id,-2,node_list,element_list,1,3,1,xpoint2)   !----------- solve for current profile (GS_inverse)
!call poisson(my_id,0,node_list,element_list,3,1,1,xpoint2)    !----------- solve for GS (from current profile)

equil_error = 0.d0
equil_value = 0.d0

if (my_id .eq. 0) then

  do i=1,node_list%n_nodes

    x_s   = node_list%node(i)%x(2,1)
    x_t   = node_list%node(i)%x(3,1)
    ps0_s = node_list%node(i)%values(1,2,1)
    ps0_t = node_list%node(i)%values(1,3,1)
    zj0_s = node_list%node(i)%values(1,2,3)
    zj0_t = node_list%node(i)%values(1,3,3)

    p_s = node_list%node(i)%values(1,2,5) * node_list%node(i)%values(1,1,6) &
        +  node_list%node(i)%values(1,1,5) * node_list%node(i)%values(1,2,6)
    p_t = node_list%node(i)%values(1,3,5) * node_list%node(i)%values(1,1,6)  &
        +  node_list%node(i)%values(1,1,5) * node_list%node(i)%values(1,3,6)

    call current(xpoint2,node_list%node(i)%x(1,1),node_list%node(i)%x(1,2), Z_xpoint, &
                 node_list%node(i)%values(1,1,1),psi_axis,psi_bnd,zjz)

!    write(*,'(A,i6,8e16.8)') ' check equil ',i,(ps0_s*zj0_t - ps0_t * zj0_s),2.d0*node_list%node(i)%x(1,1)*(x_s * p_t - x_t * p_s), &
!                                               (ps0_s*zj0_t - ps0_t * zj0_s)-2.d0*node_list%node(i)%x(1,1)*(x_s * p_t - x_t * p_s), &
!                                                zjz, node_list%node(i)%values(1,1,3)

    equil_value = equil_value + abs((ps0_s*zj0_t - ps0_t * zj0_s))
    equil_error = equil_error + abs((ps0_s*zj0_t - ps0_t * zj0_s)-2.d0*node_list%node(i)%x(1,1)*(x_s * p_t - x_t * p_s))

  enddo

  write(*,*) ' ERROR in equilibrium : ',equil_error,equil_error/equil_value

endif


if (my_id .eq. 0) then

  surface_list%n_psi =51
  if (allocated(surface_list%psi_values)) deallocate(surface_list%psi_values)
  allocate(surface_list%psi_values(surface_list%n_psi))
  do i=1,surface_list%n_psi
    surface_list%psi_values(i) =  (1. - 0.99999999 * (float(i)/float(surface_list%n_psi))**2 ) * psi_axis
  enddo

  call find_flux_surfaces(xpoint2,node_list,element_list,surface_list)! pb //
  call plot_flux_surfaces(node_list,element_list,surface_list,.true.)

  call q_profile(node_list,element_list,surface_list)

endif

return
end

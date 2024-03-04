subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use constants
use data_structure
use phys_module
use mod_poiss
use equil_info
use mod_interp, only: interp, interp_RZ
use mod_F_profile

! #define DEBUG_tauIC_initial

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm, iv, ifail, xcase2
integer    :: index0, index, n_node_start, n_index_start, j, k, ivar
real*8     :: amplitude, psi, psi_n, theta
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: R, Z, BigR, T0, BigR_s, T0_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi
real*8     :: P_ss, P_st, P_tt, R_out,Z_out,s_out,t_out 
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t,R_s, R_t, ps0_x, ps0_y, Z_s, Z_t, xjac, xjac_R, xjac_Z, direction, Btot
real*8     :: Omega, dOmega_dpsi, dOmega_dpsi2, zeta, Lam, dLam_dpsi, dLam_dpsi2
real*8     :: zn0, zT0, dn0_dpsi, dT0_dpsi, dn_dR, dn_dR2, dT_dR, dT_dR2, R2sh, rf, rf0
real*8     :: x21, x31, x41, psi2, psi3, psi4
logical    :: xpoint2
real*8     :: F_prof,          dF_dpsi,  dF_dz,  dF_dpsi2,  dF_dz2,  dF_dpsi_dz
real*8     :: FFprime_profile, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz

! --- For adding the V_ExB drift at the initial condition to balance the ion-diamagnetic drift
real*8     :: R_st, R_ss, R_tt, Z_st, Z_ss, Z_tt
real*8     :: s, t, VdiaR0, VdiaZ0, VdiaP0
real*8     :: VdiaR0_R, VdiaR0_Z, VdiaR0_RR, VdiaR0_ZZ, VdiaR0_RZ
real*8     :: VdiaZ0_R, VdiaZ0_Z, VdiaZ0_RR, VdiaZ0_ZZ, VdiaZ0_RZ
real*8     :: Vdiap0_R, Vdiap0_Z, Vdiap0_RR, Vdiap0_ZZ, Vdiap0_RZ
real*8     :: rho0, BR0, BZ0, Bp0, BB2, Ti0, BR0_R, BZ0_R, Bp0_R, BR0_Z, BZ0_Z, Bp0_Z
real*8     :: B_dot_B_R0, B_dot_B_Z0
real*8     :: AR0_p, AR0_Z, AZ0_R, AZ0_p, A30_R, A30_Z, rho0_R, rho0_Z, rho0_p, Ti0_R, Ti0_Z, Ti0_p, pi0_R, pi0_Z, pi0_p
real*8     :: A30_RZ, A30_RR, A30_ZZ, rho0_RZ, Ti0_RZ, rho0_RR, Ti0_RR, rho0_ZZ, Ti0_ZZ
real*8     :: pi0_RZ, pi0_RR, pi0_ZZ
real*8     :: F_R, F_Z
integer    :: n_s, n_t, id_left, id_right
real*8, dimension(n_var)    :: eq_g,eq_s,eq_t,eq_st,eq_ss,eq_tt

#ifdef DEBUG_tauIC_initial
real*8, dimension(n_var)    :: eq_g_check,eq_s_check,eq_t_check,eq_st_check,eq_ss_check,eq_tt_check
real*8     :: R_check, R_s_check, R_t_check, R_st_check, R_ss_check, R_tt_check, Z_check, Z_s_check, Z_t_check, Z_st_check, Z_ss_check, Z_tt_check
#endif

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (711)      *'
  write(*,*) '***************************************'
endif

! --- This old projection method should be replaced by a direct solve on the elements
! --- It is much cleaner and more generic
if ( (my_id .eq. 0) .and. (n_order .le. 3) ) then

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)

    call density(    xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                              dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature_i(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2, &
                                                              dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)

    call temperature_e(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2, &
                                                              dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)


    node_list%node(i)%values(1,:,var_uR) = 0.d0    
    node_list%node(i)%values(1,:,var_uZ) = 0.d0
    node_list%node(i)%values(1,:,var_up) = 0.d0

    node_list%node(i)%values(1,:,var_AR) = 0.d0    
    node_list%node(i)%values(1,:,var_AZ) = 0.d0

    node_list%node(i)%values(1,1,var_rho) = zn
    node_list%node(i)%values(1,2,var_rho) = dn_dpsi    * node_list%node(i)%values(1,2,var_A3)
    node_list%node(i)%values(1,3,var_rho) = dn_dpsi    * node_list%node(i)%values(1,3,var_A3)
    node_list%node(i)%values(1,4,var_rho) = dn_dpsi    * node_list%node(i)%values(1,4,var_A3)                                       &
                                          + dn_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3)


    node_list%node(i)%values(1,1,var_Ti) = zTi
    node_list%node(i)%values(1,2,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,2,var_A3)
    node_list%node(i)%values(1,3,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,3,var_A3)
    node_list%node(i)%values(1,4,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,4,var_A3)                                       &
                                        + dTi_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3) 
      
    node_list%node(i)%values(1,1,var_Te) = zTe
    node_list%node(i)%values(1,2,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,2,var_A3)
    node_list%node(i)%values(1,3,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,3,var_A3)
    node_list%node(i)%values(1,4,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,4,var_A3)                                       &
                                        + dTe_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3) 

    node_list%node(i)%deltas = 0.d0

    node_list%node(i)%psi_eq(:) = node_list%node(i)%values(1,:,1)

    ! Fprof_eq was aleady initialised in equilibrium.f90. 
    ! We fill in the values here as well, but anyway, we solve Fprof = Fprof below to ensure that the node values are clean
    ! This makes it 100% certain that all derivatives of Fprofile (when taken from the node values), will be accurate
    ! to the level of our finite elements.
    call F_profile(xpoint2, xcase2, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, &
                  F_prof,          dF_dpsi,  dF_dz,  dF_dpsi2,  dF_dz2,  dF_dpsi_dz , &
                  FFprime_profile, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz)
    node_list%node(i)%Fprof_eq(1) =   F_prof
    node_list%node(i)%Fprof_eq(2) =   dF_dpsi  * node_list%node(i)%values(1,2,var_A3) 
    node_list%node(i)%Fprof_eq(3) =   dF_dpsi  * node_list%node(i)%values(1,3,var_A3) 
    node_list%node(i)%Fprof_eq(4) = dF_dpsi    * node_list%node(i)%values(1,4,var_A3)                                             &
                                  + dF_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3) 
    if(xpoint2) then
      node_list%node(i)%values(1,2,var_rho) = node_list%node(i)%values(1,2,var_rho) + dn_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_rho) = node_list%node(i)%values(1,3,var_rho) + dn_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_rho) = node_list%node(i)%values(1,4,var_rho) + dn_dz * node_list%node(i)%x(1,4,2)           &
                                            + dn_dz2     * node_list%node(i)%x(1,2,2)            * node_list%node(i)%x(1,3,2)      &
                                            + dn_dpsi_dz * node_list%node(i)%values(1,3,var_A3)  * node_list%node(i)%x(1,2,2)      &
                                            + dn_dpsi_dz * node_list%node(i)%values(1,2,var_A3)  * node_list%node(i)%x(1,3,2)      

      node_list%node(i)%values(1,2,var_Ti) = node_list%node(i)%values(1,2,var_Ti) + dTi_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_Ti) = node_list%node(i)%values(1,3,var_Ti) + dTi_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_Ti) = node_list%node(i)%values(1,4,var_Ti) + dTi_dz * node_list%node(i)%x(1,4,2)            &
                                           + dTi_dz2     * node_list%node(i)%x(1,2,2)           * node_list%node(i)%x(1,3,2)       &
                                           + dTi_dpsi_dz * node_list%node(i)%values(1,3,var_A3) * node_list%node(i)%x(1,2,2)       &
                                           + dTi_dpsi_dz * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%x(1,3,2)      

      node_list%node(i)%values(1,2,var_Te) = node_list%node(i)%values(1,2,var_Te) + dTe_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%values(1,3,var_Te) = node_list%node(i)%values(1,3,var_Te) + dTe_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%values(1,4,var_Te) = node_list%node(i)%values(1,4,var_Te) + dTe_dz * node_list%node(i)%x(1,4,2)            &
                                           + dTe_dz2     * node_list%node(i)%x(1,2,2)           * node_list%node(i)%x(1,3,2)       &
                                           + dTe_dpsi_dz * node_list%node(i)%values(1,3,var_A3) * node_list%node(i)%x(1,2,2)       &
                                           + dTe_dpsi_dz * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%x(1,3,2)      

      node_list%node(i)%Fprof_eq(2)        = node_list%node(i)%Fprof_eq(2) + dF_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%Fprof_eq(3)        = node_list%node(i)%Fprof_eq(3) + dF_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%Fprof_eq(4)        = node_list%node(i)%Fprof_eq(4) + dF_dz * node_list%node(i)%x(1,4,2)                    &
                                           + dF_dz2     * node_list%node(i)%x(1,2,2)            * node_list%node(i)%x(1,3,2)       &
                                           + dF_dpsi_dz * node_list%node(i)%values(1,3,var_A3)  * node_list%node(i)%x(1,2,2)       &
                                           + dF_dpsi_dz * node_list%node(i)%values(1,2,var_A3)  * node_list%node(i)%x(1,3,2)      
    endif
  enddo

  if(tauIC_initial) then
#ifdef DEBUG_tauIC_initial
    open(911,file='Debug_initial_conditions.txt')
#endif
    if(n_flux .eq. 0) then
      n_s = n_radial
      n_t = n_pol
    else
      n_s = n_flux
      n_t = n_tht ! number of grids/elements in a circle
    endif
    do i_elm=1,element_list%n_elements
      ! if(i_elm > (n_s - 2) * n_t) CYCLE
      do iv=1,1 !n_vertex_max
        i = element_list%element(i_elm)%vertex(iv)
#ifdef DEBUG_tauIC_initial
        print*, "i_elm = ", i_elm, ", iv = ", iv, ", i = ", i
#endif
        psi = node_list%node(i)%values(1,1,1)
        R   = node_list%node(i)%x(1,1,1)
        Z   = node_list%node(i)%x(1,1,2)
        select case (iv)
        case(1)
          s = 0.d0; t= 0.d0
        case(2)
          s = 1.d0; t= 0.d0
        case(3)
          s = 1.d0; t= 1.d0
        case(4)
          s = 0.d0; t= 1.d0
        end select

        call interp_RZ(node_list, element_list, i_elm, s, t, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
#ifdef DEBUG_tauIC_initial
        if(iv.eq.2) then
          call interp_RZ(node_list, element_list, i_elm + n_t, 0.d0, t, R_check, R_s_check, R_t_check, R_st_check, R_ss_check, R_tt_check, Z_check, Z_s_check, Z_t_check, Z_st_check, Z_ss_check, Z_tt_check)
          if((R.ne.R_check).or.(R_s.ne.R_s_check).or.(R_t.ne.R_t_check).or.(R_st.ne.R_st_check).or.(R_ss.ne.R_ss_check).or.(R_tt.ne.R_tt_check)) then
            print*, "R = ", R, ", R_check = ", R_check
            print*, "R_s = ", R_s, ", R_s_check = ", R_s_check
            print*, "R_t = ", R_t, ", R_t_check = ", R_t_check
            print*, "R_st = ", R_st, ", R_st_check = ", R_st_check
            print*, "R_ss = ", R_ss, ", R_ss_check = ", R_ss_check
            print*, "R_tt = ", R_tt, ", R_tt_check = ", R_tt_check
          endif
          if((Z.ne.R_check).or.(Z_s.ne.Z_s_check).or.(Z_t.ne.Z_t_check).or.(Z_st.ne.Z_st_check).or.(Z_ss.ne.Z_ss_check).or.(Z_tt.ne.Z_tt_check)) then
            print*, "Z = ", Z, ", Z_check = ", Z_check
            print*, "Z_s = ", Z_s, ", Z_s_check = ", Z_s_check
            print*, "Z_t = ", Z_t, ", Z_t_check = ", Z_t_check
            print*, "Z_st = ", Z_st, ", Z_st_check = ", Z_st_check
            print*, "Z_ss = ", Z_ss, ", Z_ss_check = ", Z_ss_check
            print*, "Z_tt = ", Z_tt, ", Z_tt_check = ", Z_tt_check
          endif
        endif
#endif
        xjac = R_s * Z_t - R_t * Z_s
        if(xjac .eq. 0) xjac = 1.e19
        xjac_R  = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t   &
                + Z_st*(R_s*Z_t + R_t*Z_s)                                               &
                + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / xjac
        xjac_Z  = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s   &
                + R_st*(Z_t*R_s + Z_s*R_t)                                               &
                + Z_ss*R_t**2 - R_ss*Z_t*R_t) / xjac
        do k=1,n_var
          call interp(node_list, element_list, i_elm, k, 1, s, t, eq_g(k), eq_s(k), eq_t(k), eq_st(k), eq_ss(k), eq_tt(k))
#ifdef DEBUG_tauIC_initial
          if(iv.eq.2) then
            call interp(node_list, element_list, i_elm + n_t, k, 1, 0.d0, t, eq_g_check(k), eq_s_check(k), eq_t_check(k), eq_st_check(k), eq_ss_check(k), eq_tt_check(k))
            if((eq_g(k).ne.eq_g_check(k)).or.(eq_s(k).ne.eq_s_check(k)).or.(eq_t(k).ne.eq_t_check(k)).or.(eq_st(k).ne.eq_st_check(k)).or.(eq_ss(k).ne.eq_ss_check(k)).or.(eq_tt(k).ne.eq_tt_check(k))) then
              print*, "k = ", k
              print*, "eq_g(k) = ", eq_g(k), ", eq_g_check(k) = ", eq_g_check(k)
              print*, "eq_s(k) = ", eq_s(k), ", eq_s_check(k) = ", eq_s_check(k)
              print*, "eq_t(k) = ", eq_t(k), ", eq_t_check(k) = ", eq_t_check(k)
              print*, "eq_st(k) = ", eq_st(k), ", eq_st_check(k) = ", eq_st_check(k)
              print*, "eq_ss(k) = ", eq_ss(k), ", eq_ss_check(k) = ", eq_ss_check(k)
              print*, "eq_tt(k) = ", eq_tt(k), ", eq_tt_check(k) = ", eq_tt_check(k)
            endif
          endif
#endif
        enddo
#ifdef DEBUG_tauIC_initial
        if((R.ne.node_list%node(i)%x(1,1,1)).or.((Z.ne.node_list%node(i)%x(1,1,2))).or.(xjac.eq.0)) then ! .or.(mod(i,1000) <= 4)
          write(911,*) "i_elm = ", i_elm, ", iv = ", iv, ", i = ", i, ", ERROR!!!"
          write(911,*) "xjac=", xjac

          write(911,*) "R   =", R  , ", x(1,1)=", node_list%node(i)%x(1,1,1)
          write(911,*) "R_s =", R_s, ", x(2,1)=", node_list%node(i)%x(1,2,1)
          write(911,*) "R_t =", R_t, ", x(3,1)=", node_list%node(i)%x(1,3,1)
          write(911,*) "R_st=", R_st,", x(4,1)=", node_list%node(i)%x(1,4,1)
          
          write(911,*) "Z   =", Z  , ", x(1,2)=", node_list%node(i)%x(1,1,2)
          write(911,*) "Z_s =", Z_s, ", x(2,2)=", node_list%node(i)%x(1,2,2)
          write(911,*) "Z_t =", Z_t, ", x(3,2)=", node_list%node(i)%x(1,3,2)
          write(911,*) "Z_st=", Z_st,", x(4,2)=", node_list%node(i)%x(1,4,2)
          
          write(911,*) "rho   =", eq_g(var_rho)  , ", rho(1)=", node_list%node(i)%values(1,1,var_rho)
          write(911,*) "rho_s =", eq_s(var_rho)  , ", rho(2)=", node_list%node(i)%values(1,2,var_rho)
          write(911,*) "rho_t =", eq_t(var_rho)  , ", rho(3)=", node_list%node(i)%values(1,3,var_rho)
          write(911,*) "rho_st=", eq_st(var_rho) , ", rho(4)=", node_list%node(i)%values(1,4,var_rho)

          write(911,*) "psi   =", eq_g(var_A3)  , ", A3(1)=", node_list%node(i)%values(1,1,var_A3)
          write(911,*) "psi_s =", eq_s(var_A3)  , ", A3(2)=", node_list%node(i)%values(1,2,var_A3)
          write(911,*) "psi_t =", eq_t(var_A3)  , ", A3(3)=", node_list%node(i)%values(1,3,var_A3)
          write(911,*) "psi_st=", eq_st(var_A3) , ", A3(4)=", node_list%node(i)%values(1,4,var_A3)
        endif
#endif
        AR0_p = 0.d0
        AZ0_p = 0.d0
        AZ0_R = (  Z_t * eq_s(var_AZ) - Z_s * eq_t(var_AZ)) / xjac ! = 0.d0
        AR0_Z = (- R_t * eq_s(var_AR) + R_s * eq_t(var_AR)) / xjac ! = 0.d0
        A30_R = (  Z_t * eq_s(var_A3) - Z_s * eq_t(var_A3)) / xjac
        A30_Z = (- R_t * eq_s(var_A3) + R_s * eq_t(var_A3)) / xjac

        A30_RZ = (- eq_ss(var_A3) * Z_t*R_t - eq_tt(var_A3) * R_s*Z_s                     &
              + eq_st(var_A3) * (Z_s*R_t  + Z_t*R_s  ) &
              - eq_s(var_A3)  * (R_st*Z_t - R_tt*Z_s ) &
              - eq_t(var_A3) * (R_st*Z_s  - R_ss*Z_t ) ) / xjac**2             &
              - xjac_R * (- eq_s(var_A3) * R_t + eq_t(var_A3) * R_s )   / xjac**2

        A30_RR = (eq_ss(var_A3) * Z_t**2 - 2.d0*eq_st(var_A3) * Z_s*Z_t + eq_tt(var_A3) * Z_s**2 &
              + eq_s(var_A3) * (Z_st*Z_t - Z_tt*Z_s ) &
              + eq_t(var_A3) * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac**2             &
              - xjac_R * (eq_s(var_A3)* Z_t - eq_t(var_A3) * Z_s)  / xjac**2

        A30_ZZ = (eq_ss(var_A3) * R_t**2 - 2.d0*eq_st(var_A3) * R_s*R_t + eq_tt(var_A3) * R_s**2 &
              + eq_s(var_A3) * (R_st*R_t - R_tt*R_s ) &
              + eq_t(var_A3) * (R_st*R_s - R_ss*R_t ) ) / xjac**2             &
              - xjac_Z * (- eq_s(var_A3) * R_t + eq_t(var_A3) * R_s )  / xjac**2

        rho0   = max(eq_g(var_rho), 1d-9)
        rho0_R = (  Z_t * eq_s(var_rho) - Z_s * eq_t(var_rho)) / xjac
        rho0_Z = (- R_t * eq_s(var_rho) + R_s * eq_t(var_rho)) / xjac
        rho0_p = 0.d0
        rho0_RZ = (- eq_ss(var_rho) * Z_t*R_t - eq_tt(var_rho) * R_s*Z_s                     &
              + eq_st(var_rho) * (Z_s*R_t  + Z_t*R_s  ) &
              - eq_s(var_rho)  * (R_st*Z_t - R_tt*Z_s ) &
              - eq_t(var_rho) * (R_st*Z_s  - R_ss*Z_t ) ) / xjac**2             &
              - xjac_R * (- eq_s(var_rho) * R_t + eq_t(var_rho) * R_s )   / xjac**2

        rho0_RR = (eq_ss(var_rho) * Z_t**2 - 2.d0*eq_st(var_rho) * Z_s*Z_t + eq_tt(var_rho) * Z_s**2 &
              + eq_s(var_rho) * (Z_st*Z_t - Z_tt*Z_s ) &
              + eq_t(var_rho) * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac**2             &
              - xjac_R * (eq_s(var_rho)* Z_t - eq_t(var_rho) * Z_s)  / xjac**2
              
        rho0_ZZ = (eq_ss(var_rho) * R_t**2 - 2.d0*eq_st(var_rho) * R_s*R_t + eq_tt(var_rho) * R_s**2 &
              + eq_s(var_rho) * (R_st*R_t - R_tt*R_s ) &
              + eq_t(var_rho) * (R_st*R_s - R_ss*R_t ) ) / xjac**2             &
              - xjac_Z * (- eq_s(var_rho) * R_t + eq_t(var_rho) * R_s )  / xjac**2

        Ti0   = eq_g(var_Ti)
        Ti0_R = (  Z_t * eq_s(var_Ti) - Z_s * eq_t(var_Ti)) / xjac
        Ti0_Z = (- R_t * eq_s(var_Ti) + R_s * eq_t(var_Ti)) / xjac
        Ti0_p = 0.d0
        Ti0_RZ = (- eq_ss(var_Ti) * Z_t*R_t - eq_tt(var_Ti) * R_s*Z_s                     &
              + eq_st(var_Ti) * (Z_s*R_t  + Z_t*R_s  ) &
              - eq_s(var_Ti)  * (R_st*Z_t - R_tt*Z_s ) &
              - eq_t(var_Ti) * (R_st*Z_s  - R_ss*Z_t ) ) / xjac**2             &
              - xjac_R * (- eq_s(var_Ti) * R_t + eq_t(var_Ti) * R_s )   / xjac**2
        
        Ti0_RR = (eq_ss(var_Ti) * Z_t**2 - 2.d0*eq_st(var_Ti) * Z_s*Z_t + eq_tt(var_Ti) * Z_s**2 &
              + eq_s(var_Ti) * (Z_st*Z_t - Z_tt*Z_s ) &
              + eq_t(var_Ti) * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac**2             &
              - xjac_R * (eq_s(var_Ti)* Z_t - eq_t(var_Ti) * Z_s)  / xjac**2

        Ti0_ZZ = (eq_ss(var_Ti) * R_t**2 - 2.d0*eq_st(var_Ti) * R_s*R_t + eq_tt(var_Ti) * R_s**2 &
              + eq_s(var_Ti) * (R_st*R_t - R_tt*R_s ) &
              + eq_t(var_Ti) * (R_st*R_s - R_ss*R_t ) ) / xjac**2             &
              - xjac_Z * (- eq_s(var_Ti) * R_t + eq_t(var_Ti) * R_s )  / xjac**2

        BR0 = ( A30_Z - AZ0_p ) / R
        BZ0 = ( AR0_p - A30_R ) / R
        Bp0 = ( AZ0_R - AR0_Z ) + F_prof / R

        F_R = dF_dpsi * A30_R
        F_Z = dF_dpsi * A30_Z

        BR0_R = -A30_Z / R**2 + A30_RZ / R
        BZ0_R =  A30_R / R**2 - A30_RR / R
        Bp0_R = -F_prof / R**2 + F_R / R

        BR0_Z = A30_ZZ / R
        BZ0_Z = -A30_RZ / R
        Bp0_Z = F_Z / R

        B_dot_B_R0 = BR0 * BR0_R + BZ0 * BZ0_R + Bp0 * Bp0_R
        B_dot_B_Z0 = BR0 * BR0_Z + BZ0 * BZ0_Z + Bp0 * Bp0_Z

        BB2   = Bp0**2 + BR0**2 + BZ0**2
        pi0_R = rho0_R * Ti0 + rho0 * Ti0_R
        pi0_Z = rho0_Z * Ti0 + rho0 * Ti0_Z
        pi0_p = rho0_p * Ti0 + rho0 * Ti0_p

        VdiaR0 = - tauIC * F0 / (R * rho0 * BB2) * (  BZ0*pi0_p - R*Bp0*pi0_Z)
        VdiaZ0 = - tauIC * F0 / (R * rho0 * BB2) * (R*BP0*pi0_R -   BR0*pi0_p)
        VdiaP0 = - tauIC * F0 / (    rho0 * BB2) * (  BR0*pi0_Z -   BZ0*pi0_R)

        pi0_RZ = rho0_R * Ti0_Z + rho0 * Ti0_RZ + Ti0_R * rho0_Z + Ti0 * rho0_RZ
        pi0_ZZ = rho0_Z * Ti0_Z + rho0 * Ti0_ZZ + Ti0_Z * rho0_Z + Ti0 * rho0_ZZ
        pi0_RR = rho0_R * Ti0_R + rho0 * Ti0_RR + Ti0_R * rho0_R + Ti0 * rho0_RR

        VdiaR0_R = tauIC * F0 * ((Ti0 * rho0_Z + rho0 * Ti0_Z) * (F_R / (R * rho0 * BB2) - F_prof / (R * rho0 * BB2)**2 * (rho0 * BB2 + R * BB2 * rho0_R + 2 * R * rho0 * B_dot_B_R0)) &
                + F_prof / (R * rho0 * BB2) * (Ti0_R * rho0_Z + rho0_R * Ti0_Z + Ti0 * rho0_RZ + rho0 * Ti0_RZ)) 
        VdiaR0_Z = tauIC * F0 * ((Ti0 * rho0_Z + rho0 * Ti0_Z) * (F_Z / (R * rho0 * BB2) - F_prof / (R * rho0 * BB2)**2 * (             R * BB2 * rho0_Z + 2 * R * rho0 * B_dot_B_Z0)) &
                + F_prof / (R * rho0 * BB2) * (2 * Ti0_Z * rho0_Z              + Ti0 * rho0_ZZ + rho0 * Ti0_ZZ))
        
        vdiaR0_RR = 0.d0
        vdiaR0_RZ = 0.d0
        vdiaR0_ZZ = 0.d0
        
        VdiaZ0_R =-tauIC * F0 * ((Ti0 * rho0_R + rho0 * Ti0_R) * (F_R / (R * rho0 * BB2) - F_prof / (R * rho0 * BB2)**2 * (rho0 * BB2 + R * BB2 * rho0_R + 2 * R * rho0 * B_dot_B_R0)) &
                + F_prof / (R * rho0 * BB2) * (2 * Ti0_R * rho0_R              + Ti0 * rho0_RR + rho0 * Ti0_RR))
        VdiaZ0_Z =-tauIC * F0 * ((Ti0 * rho0_R + rho0 * Ti0_R) * (F_Z / (R * rho0 * BB2) - F_prof / (R * rho0 * BB2)**2 * (             R * BB2 * rho0_Z + 2 * R * rho0 * B_dot_B_Z0)) &
                + F_prof / (R * rho0 * BB2) * (Ti0_Z * rho0_R + rho0_Z * Ti0_R + Ti0 * rho0_RZ + rho0 * Ti0_RZ))
        
        vdiaZ0_RR = 0.d0
        vdiaZ0_RZ = 0.d0
        vdiaZ0_ZZ = 0.d0

        Vdiap0_R = -tauIC * F0 * ((A30_Z * pi0_Z + A30_R * pi0_R) * -1.d0 / (R * rho0 * BB2)**2 * (rho0 * BB2 + R * BB2 * rho0_R + 2 * R * rho0 * B_dot_B_R0) &
                  + 1.d0  / (R * rho0 * BB2) * (A30_RZ * pi0_Z + A30_Z * pi0_RZ + A30_RR * pi0_R + A30_R * pi0_RR))
        Vdiap0_Z = -tauIC * F0 * ((A30_Z * pi0_Z + A30_R * pi0_R) * -1.d0 / (R * rho0 * BB2)**2 * (             R * BB2 * rho0_Z + 2 * R * rho0 * B_dot_B_Z0) &
                  + 1.d0  / (R * rho0 * BB2) * (A30_ZZ * pi0_Z + A30_Z * pi0_ZZ + A30_RZ * pi0_R + A30_R * pi0_RZ))

        vdiap0_RR = 0.d0
        vdiap0_RZ = 0.d0
        vdiap0_ZZ = 0.d0
        
        node_list%node(i)%values(1,1,var_uR) = VdiaR0
        node_list%node(i)%values(1,2,var_uR) = VdiaR0_R * R_s + VdiaR0_Z * Z_s
        node_list%node(i)%values(1,3,var_uR) = VdiaR0_R * R_t + VdiaR0_Z * Z_t
        ! node_list%node(i)%values(1,4,var_uR) = VdiaR0_R * R_st + VdiaR0_Z * Z_st   &
        !                                      + vdiaR0_RR * R_s * R_t               &
        !                                      + VdiaR0_RZ * (R_s * Z_t + Z_s * R_t) &
        !                                      + VdiaR0_ZZ * Z_s * Z_t

        node_list%node(i)%values(1,1,var_uZ) = VdiaZ0
        node_list%node(i)%values(1,2,var_uZ) = VdiaZ0_R * R_s + VdiaZ0_Z * Z_s
        node_list%node(i)%values(1,3,var_uZ) = VdiaZ0_R * R_t + VdiaZ0_Z * Z_t
        ! node_list%node(i)%values(1,4,var_uZ) = VdiaZ0_R * R_st + VdiaZ0_Z * Z_st   &
        !                                      + vdiaZ0_RR * R_s * R_t               &
        !                                      + VdiaZ0_RZ * (R_s * Z_t + Z_s * R_t) &
        !                                      + VdiaZ0_ZZ * Z_s * Z_t
        
        node_list%node(i)%values(1,1,var_up) = VdiaP0
        node_list%node(i)%values(1,2,var_up) = Vdiap0_R * R_s + Vdiap0_Z * Z_s
        node_list%node(i)%values(1,3,var_up) = Vdiap0_R * R_t + Vdiap0_Z * Z_t
        ! node_list%node(i)%values(1,4,var_up) = Vdiap0_R * R_st + Vdiap0_Z * Z_st   &
        !                                      + vdiap0_RR * R_s * R_t               &
        !                                      + Vdiap0_RZ * (R_s * Z_t + Z_s * R_t) &
        !                                      + Vdiap0_ZZ * Z_s * Z_t

        if(xpoint2) then
          write(*,*) 'Warnning, the partial_z(2) derivatives for ion-diamagnetic drift equlibria are not implemented yet, please check!!!!!!'
        endif
      enddo
    enddo
    ! calculate the d2/dsdt using finite difference method
    do i=1,node_list%n_nodes
      if(i <= n_t .or. i > ((n_s - 1) * n_t)) CYCLE ! skip the innermost and outermost circles
      id_right = i + 1
      id_left  = i - 1
      if(mod(i,n_t) .eq. 1) then
        id_left = id_left + n_t
      elseif(mod(i,n_t) .eq. 0) then
        id_right = id_right - n_t
      endif
      node_list%node(i)%values(1,4,var_uR) = ((node_list%node(id_right)%values(1,2,var_uR) - node_list%node(id_left)%values(1,2,var_uR))/2.d0 &
                                            + (node_list%node(i + n_t)%values(1,3,var_uR) - node_list%node(i - n_t)%values(1,3,var_uR))/2.d0)/2.d0

      node_list%node(i)%values(1,4,var_uZ) = ((node_list%node(id_right)%values(1,2,var_uZ) - node_list%node(id_left)%values(1,2,var_uZ))/2.d0 &
                                            + (node_list%node(i + n_t)%values(1,3,var_uZ) - node_list%node(i - n_t)%values(1,3,var_uZ))/2.d0)/2.d0

      node_list%node(i)%values(1,4,var_up) = ((node_list%node(id_right)%values(1,2,var_up) - node_list%node(id_left)%values(1,2,var_up))/2.d0 &
                                            + (node_list%node(i + n_t)%values(1,3,var_up) - node_list%node(i - n_t)%values(1,3,var_up))/2.d0)/2.d0
    enddo
#ifdef DEBUG_tauIC_initial
    close(911)
#endif
  endif ! if tauIC_initial
endif ! if n_order<=3

! --- Variable projection is better at higher order...
! --- (by the way, we could use this for n_order=3 and remove all the above as well, 
! --- and remove all derivatives from profiles functions, which are not really needed, 
! --- except dn_dpsi and dT_dpsi for current profile...)
if (n_order .ge. 5) then
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_A3,var_rho,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_A3,var_Ti,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_A3,var_Te,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
endif



! --- This is the special Poisson for Fprofile (it will not overwrite var_A3)
call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
             var_A3,710,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)      ! inverse Poisson


!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,1)
      Z   = node_list%node(i)%x(1,1,2)
      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

     
      ! Initialise perturbation for A3 nonzero n harmonics
      node_list%node(i)%values(in,:,:)= 0.d0

      node_list%node(i)%values(in,1,var_A3) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,var_A3) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,2,var_A3)
      node_list%node(i)%values(in,3,var_A3) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,3,var_A3)
      node_list%node(i)%values(in,4,var_A3) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,4,var_A3)

      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. UPPER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. LOWER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

enddo




return
end

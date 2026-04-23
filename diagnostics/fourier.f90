!> Allows to determine the poloidal magnetic coordinate theta_mag from field line tracing
!! in the n=0 component of the magnetic field and to perform a Fourier transformation of
!! physical quantities in (theta_mag, phi)
module fourier
  
  use constants
  use tr_module 
  use mod_parameters,      only: n_vertex_max, n_degrees, n_plane, n_tor, n_var, variable_names
  use nodes_elements,  only: node_list, element_list
  use phys_module,     only: F0, xpoint, xcase
  use mod_straight_field_line

  implicit none
  
  save
  
  public
  
  integer, parameter   :: FFTW_ESTIMATE = 64  !< (constant of the FFTW library)


  
  contains
  
  
  
  !> Fourier-transform the physical variables in the magnetic angle theta_mag.
  subroutine transform_qttys(mapping, vfour, nTht)
    use mod_basisfunctions
#ifdef fullmhd
    use mod_model_settings, only: var_A3, var_AR, var_AZ, var_uR, var_uZ, var_up, var_rho, var_T, var_Ti, var_Te, with_TiTe
    use mod_math_operators
    use mod_interp
#endif
    use basis_at_gaussian
    type(t_theta_mapping), intent(in)    :: mapping
    complex, allocatable,  intent(inout) :: vfour(:,:,:,:) !< Transformed quantities (m,n,irad,ivar)
    integer,               intent(in)    :: nTht !< Number of poloidal equidistant points
    
    real*8, allocatable :: vve(:,:,:,:) ! Variable values (ipol,itor,irad,ivar)
    real*8, dimension(4,n_degrees)  :: G, G_s, G_t, G_st, G_ss, G_tt
    integer :: i, j, k, l, iharm, nn, kv, iv, kf
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail
    real*8  :: v, basis_function
    integer*8 :: fftw_plan
    integer :: nequidist_tor, nequidist_pts
    
#ifdef fullmhd
    real*8 :: R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
    real*8 :: Ps0,Ps0_s,Ps0_t,Ps0_st,Ps0_ss,Ps0_tt
    real*8 :: eq,eq_s,eq_t,eq_st,eq_ss,eq_tt
    real*8 :: F_prof, F_prof_R, F_prof_Z
    real*8 :: ps0_x, ps0_y, psi_abs, wR2r, wZ2r, wR2p, wZ2p
    real*8 :: AR_p, AZ_p, A3_R, A3_Z, AR_Z, AZ_R
    real*8 :: A3_RR, A3_RZ, A3_ZZ, A3_Rp, A3_Zp, A3_p, A3_pp
    real*8 :: AR_RR, AR_RZ, AR_Rp, AR_Zp, AR_pp, AR_ZZ
    real*8 :: AZ_RR, AZ_RZ, AZ_Rp, AZ_Zp, AZ_pp, AZ_ZZ
    real*8 :: rho, rho_R, rho_Z, rho_p
    real*8 :: T, T_s, T_t, T_R, T_Z
    real*8 :: UR_R, UZ_Z, Up_p
    real*8 :: A_radial, A_poloidal, U_radial, U_poloidal, B_radial, B_poloidal
    real*8 :: BR, BZ, Bp
    real*8 :: BR_R, BR_Z, BR_p, BZ_R, BZ_Z, BZ_p
    real*8 :: BP_R, BP_Z, BP_p
    real*8 :: JR, JZ, Jphi
    real*8 :: gradP_R, gradP_Z
    real*8 :: JxB_radial, gradP_radial, RUxB_phi
    real*8 :: rho_divU, U_gradRho
    real*8 :: xjac, xjac_R, xjac_Z
#endif

    nequidist_tor = n_plane
    nequidist_pts =  nTht
    
    if ( allocated(vfour) ) deallocate(vfour)
#ifdef fullmhd
#define n_var_additional 14
    call initialise_basis
    allocate(vve(nequidist_pts,nequidist_tor,mapping%nstpts,n_var + n_var_additional))
    allocate(vfour(nequidist_pts/2+1,nequidist_tor,mapping%nstpts,n_var + n_var_additional))
#else
    allocate(vve(nequidist_pts,nequidist_tor,mapping%nstpts,n_var))
    allocate(vfour(nequidist_pts/2+1,nequidist_tor,mapping%nstpts,n_var))
#endif
    vve   = 0.d0
    vfour = 0.d0
    
    do k = 1, mapping%nstpts ! radial positions
      
      write(*,'(1x,a,i4)') 'Transforming variables on surface', k
      
      do j = 1, nequidist_tor  ! toroidal positions
#ifdef fullmhd
        !$omp parallel do                                                                          &
        !$omp   default(shared)                                                                    &
        !$omp   firstprivate(iharm,R_out,Z_out,i_elm_out,s_out,t_out,ifail,v,l,nn,kv,iv,kf,        &
        !$omp     basis_function,G,G_s,G_t,G_st,G_ss,G_tt,                                         &
        !$omp     R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt,                               &
        !$omp     Ps0,Ps0_s,Ps0_t,Ps0_st,Ps0_ss,Ps0_tt,                                            &
        !$omp     eq,eq_s,eq_t,eq_st,eq_ss,eq_tt,                                                  &
        !$omp     F_prof, F_prof_R, F_prof_Z,                                                      &
        !$omp     A3_R, A3_Z, A3_p, A3_RR, A3_RZ, A3_ZZ, A3_Rp, A3_Zp, A3_pp,                      &
        !$omp     AR_p, AR_Z, AR_pp, AR_Rp, AR_Zp, AR_RR, AR_RZ, AR_ZZ,                            &
        !$omp     AZ_p, AZ_R, AZ_pp, AZ_Rp, AZ_Zp, AZ_RR, AZ_RZ, AZ_ZZ,                            &
        !$omp     BR_R, BR_Z, BR_p, BZ_R, BZ_Z, BZ_p, BP_R, BP_Z, BP_p,                            &
        !$omp     T, T_s, T_t, T_R, T_Z,                                                           &
        !$omp     rho, rho_R, rho_Z, rho_p,                                                        &
        !$omp     UR_R, UZ_Z, Up_p,                                                                &
        !$omp     JR, JZ, Jphi,                                                                    &
        !$omp     gradP_R, gradP_Z, JxB_radial, gradP_radial, RUxB_phi,                            &
        !$omp     rho_divU, U_gradRho,                                                             &
        !$omp     ps0_x, ps0_y, psi_abs, wR2r, wZ2r, wR2p, wZ2p,                                   &
        !$omp     A_radial, A_poloidal, U_radial, U_poloidal, B_radial, B_poloidal,                &
        !$omp     BR, BZ, Bp,                                                                      &
        !$omp     xjac, xjac_R, xjac_Z)                                                            &
        !$omp   private(i)
#else
        !$omp parallel do                                                                          &
        !$omp   default(shared)                                                                    &
        !$omp   firstprivate(iharm,R_out,Z_out,i_elm_out,s_out,t_out,ifail,v,l,nn,kv,iv,kf,        &
        !$omp     basis_function,G,G_s,G_t,G_st,G_ss,G_tt)                                         &
        !$omp   private(i)
#endif
        do i = 1, nequidist_pts  ! poloidal positions
          
          call find_RZ(node_list,element_list,mapping%rre(k,i-1),mapping%zze(k,i-1),R_out,Z_out,   &
            i_elm_out,s_out,t_out,ifail)
#ifdef fullmhd
          wR2r = 0.d0; wZ2r = 0.d0; wR2p = 0.d0; wZ2p = 0.d0
          AR_p = 0.d0; AZ_p = 0.d0; A3_R = 0.d0; A3_Z = 0.d0; AR_Z = 0.d0; AZ_R = 0.d0
          A3_RR = 0.d0; A3_RZ = 0.d0; A3_ZZ = 0.d0; A3_Rp = 0.d0; A3_Zp = 0.d0
          A3_p = 0.d0; A3_pp = 0.d0
          AR_RR = 0.d0; AR_RZ = 0.d0; AR_Rp = 0.d0; AR_Zp = 0.d0; AR_pp = 0.d0; AR_ZZ = 0.d0
          AZ_RR = 0.d0; AZ_RZ = 0.d0; AZ_Rp = 0.d0; AZ_Zp = 0.d0; AZ_pp = 0.d0; AZ_ZZ = 0.d0
          F_prof_R = 0.d0; F_prof_Z = 0.d0
          rho = 0.d0; rho_R = 0.d0; rho_Z = 0.d0; rho_p = 0.d0
          T = 0.d0; T_R = 0.d0; T_Z = 0.d0
          UR_R = 0.d0; UZ_Z = 0.d0; Up_p = 0.d0
          call interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
          xjac = op_xjac(R_s, R_t, Z_s, Z_t)
          if(xjac .eq. 0.d0) then
            xjac = 1.d-12
          end if
          xjac_R = op_xjac_R(R_s, R_t, R_st, R_ss, R_tt,              &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac)
          xjac_Z = op_xjac_Z(R_s, R_t, R_st, R_ss, R_tt,              &
                        Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                        xjac)
          call interp(node_list,element_list,i_elm_out,var_A3,1,s_out,t_out,Ps0,Ps0_s,Ps0_t,Ps0_st,Ps0_ss,Ps0_tt)
          if (xjac .gt. 1.d-6) then
            ps0_x = op_grad_R(Ps0_s, Ps0_t, R_s, R_t, Z_s, Z_t, xjac)
            ps0_y = op_grad_Z(Ps0_s, Ps0_t, R_s, R_t, Z_s, Z_t, xjac)
            psi_abs = dsqrt(ps0_x * ps0_x + ps0_y * ps0_y)
            wR2r =   ps0_x / psi_abs
            wZ2r =   ps0_y / psi_abs
            wR2p = - ps0_y / psi_abs
            wZ2p =   ps0_x / psi_abs
          endif
#endif
          
          do iharm = 1, n_tor    ! toroidal harmonics

            nn = iharm / 2 ! toroidal mode number (without periodicity)

            if ( iharm == 1 ) then
	            basis_function = 1.d0
            else if ( MOD(iharm,2) == 0 ) then
	            basis_function = cos(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            else
	            basis_function = sin(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            end if
            if(abs(basis_function - HZ(iharm,j)) .gt. 1.d-14) then
              write(*,*) 'Warning: basis function from direct calc and from HZ differ at iharm=', iharm, ' toroidal position', j, 'by', abs(basis_function - HZ(iharm,j))
              write(*,*) 'This probably is caused by nequidist_tor not equal to n_plane.'
              write(*,*) 'In such case, the basis function defined above is inconsistent with HZ_p, HZ_pp used below.'
            end if
            if(nequidist_tor .eq. n_plane) then
              basis_function = HZ(iharm,j)
            endif
#ifdef fullmhd
            do l = 1, n_var
              call interp(node_list,element_list,i_elm_out,l,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
              !$omp atomic
              vve(i,j,k,l) = vve(i,j,k,l) + eq * basis_function
            enddo
#else	    
            call basisfunctions(s_out,t_out,G,G_s,G_t,G_st,G_ss,G_tt)
  
	          do l = 1, n_var

              v = 0.d0
  
              do kv = 1, n_vertex_max  ! 4 vertices
    
                iv = element_list%element(i_elm_out)%vertex(kv)  ! the node number
    
                do kf = 1, n_degrees       ! basis functions
    
                  v = v + node_list%node(iv)%values(iharm,kf,l)                                    &
                    * element_list%element(i_elm_out)%size(kv,kf) * G(kv,kf)
    
                end do
    
              end do

              !$omp atomic
              vve(i,j,k,l) = vve(i,j,k,l) + v * basis_function
	          end do
#endif
#ifdef fullmhd
            call interp(node_list,element_list,i_elm_out,var_A3,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            A3_Z = A3_Z + op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            A3_R = A3_R + op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            A3_p  = A3_p  + eq * HZ_p(iharm,j)
            A3_pp = A3_pp + eq * HZ_pp(iharm,j)

            A3_RZ = A3_RZ + op_lap_RZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function
            A3_ZZ = A3_ZZ + op_lap_ZZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function
            A3_Zp = A3_Zp + ( - R_t * eq_s + R_s * eq_t ) / xjac * HZ_p(iharm,j)
            A3_RR = A3_RR + op_lap_RR(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function
            A3_Rp = A3_Rp + (   Z_t * eq_s - Z_s * eq_t ) / xjac * HZ_p(iharm,j)

            call interp(node_list,element_list,i_elm_out,var_AR,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            AR_p = AR_p + eq * HZ_p(iharm,j)
            AR_Z = AR_Z + op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            AR_pp = AR_pp + eq * HZ_pp(iharm,j)

            AR_Rp = AR_Rp + (   Z_t * eq_s - Z_s * eq_t ) / xjac * HZ_p(iharm,j)
            AR_Zp = AR_Zp + ( - R_t * eq_s + R_s * eq_t ) / xjac * HZ_p(iharm,j)
            AR_RZ = AR_RZ + op_lap_RZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function
            AR_ZZ = AR_ZZ + op_lap_ZZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function

            call interp(node_list,element_list,i_elm_out,var_AZ,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            AZ_p = AZ_p + eq * HZ_p(iharm,j)
            AZ_R = AZ_R + op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            AZ_pp = AZ_pp + eq * HZ_pp(iharm,j)
            AZ_Rp = AZ_Rp + (   Z_t * eq_s - Z_s * eq_t ) / xjac * HZ_p(iharm,j)
            AZ_Zp = AZ_Zp + ( - R_t * eq_s + R_s * eq_t ) / xjac * HZ_p(iharm,j)
            AZ_RR = AZ_RR + op_lap_RR(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function
            AZ_RZ = AZ_RZ + op_lap_RZ(eq_s, eq_t, eq_st, eq_ss, eq_tt,              &
                                      R_s, R_t, R_st, R_ss, R_tt,                   &
                                      Z_s, Z_t, Z_st, Z_ss, Z_tt,                   &
                                      xjac, xjac_R, xjac_Z) * basis_function

            call interp(node_list,element_list,i_elm_out,var_rho,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            rho   = rho + eq * basis_function
            rho_R = rho_R + op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            rho_Z = rho_Z + op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            rho_p = rho_p + eq * HZ_p(iharm,j)
            if (with_TiTe) then
              call interp(node_list,element_list,i_elm_out,var_Ti,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
              T   = T + eq * basis_function
              T_s = eq_s
              T_t = eq_t
              call interp(node_list,element_list,i_elm_out,var_Te,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
              T   = T + eq * basis_function
              T_s = T_s + eq_s
              T_t = T_t + eq_t
            else
              call interp(node_list,element_list,i_elm_out,var_T,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
              T   = T + eq * basis_function
              T_s = eq_s
              T_t = eq_t
            endif
            T_R = T_R + op_grad_R(T_s, T_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            T_Z = T_Z + op_grad_Z(T_s, T_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function

            call interp(node_list,element_list,i_elm_out,var_uR,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            UR_R = UR_R + op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function
            
            call interp(node_list,element_list,i_elm_out,var_uZ,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            UZ_Z = UZ_Z + op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * basis_function

            call interp(node_list,element_list,i_elm_out,var_up,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            Up_p = Up_p + eq * HZ_p(iharm,j)
#endif            
          end do ! iharm
#ifdef fullmhd
          call interp(node_list,element_list,i_elm_out,710,1,s_out,t_out,F_prof,eq_s,eq_t,eq_st,eq_ss,eq_tt)
          F_prof_R      = (   Z_t * eq_s   - Z_s * eq_t   ) / xjac
          F_prof_Z      = ( - R_t * eq_s   + R_s * eq_t   ) / xjac
          BR = ( A3_Z - AZ_p )/ R
          BZ = ( AR_p - A3_R )/ R
          BP = ( AZ_R - AR_Z )    + F_prof/ R
          BR_R = -1.d0 / (R * R) * (A3_Z - AZ_p) + (A3_RZ - AZ_Rp) / R
          BR_Z = (A3_ZZ - AZ_Zp) / R
          BR_p = (A3_Zp - AZ_pp) / R
          BZ_R = -1.d0 / (R * R) * (AR_p - A3_R) + (AR_Rp - A3_RR) / R
          BZ_Z = (AR_Zp - A3_RZ) / R
          BZ_p = (AR_pp - A3_Rp) / R
          BP_R = (AZ_RR - AR_RZ) + F_prof_R / R - F_prof / (R * R)
          BP_Z = (AZ_RZ - AR_ZZ) + F_prof_Z / R
          BP_p = (AZ_Rp - AR_Zp)
          JR = (R * BP_Z - BZ_p) / R
          JZ = (BR_p - R * BP_R - BP) / R
          Jphi = BZ_R - BR_Z

          GradP_R    = rho_R*T + rho*T_R
          GradP_Z    = rho_Z*T + rho*T_Z

          if(abs(rho - vve(i,j,k,var_rho)) .gt. 1.d-14) then
            write(*,*) 'Warning: density from interpolation is different from vve at iharm=', iharm, ' toroidal position', j, 'and poloidal position', i, 'by', abs(rho - vve(i,j,k,var_rho))
          endif
          
          JxB_radial = ((JZ * BP - Jphi * BZ) * wR2r + (Jphi * BR - JR * BP) * wZ2r)
          gradP_radial = (gradP_R * wR2r + gradP_Z * wZ2r)
          RUxB_phi = (vve(i,j,k,var_UR) * BZ - vve(i,j,k,var_UZ) * BR) * R
          rho_divU = rho * (vve(i,j,k,var_UR) / R + UR_R + UZ_Z + Up_p / R)
          U_gradRho = vve(i,j,k,var_UR) * rho_R + vve(i,j,k,var_UZ) * rho_Z + vve(i,j,k,var_up) * rho_p / R

          A_radial   = vve(i,j,k,var_AR) * wR2r + vve(i,j,k,var_AZ) * wZ2r
          A_poloidal = vve(i,j,k,var_AR) * wR2p + vve(i,j,k,var_AZ) * wZ2p
          U_radial   = vve(i,j,k,var_UR) * wR2r + vve(i,j,k,var_UZ) * wZ2r
          U_poloidal = vve(i,j,k,var_UR) * wR2p + vve(i,j,k,var_UZ) * wZ2p
          B_radial   = BR                * wR2r + BZ                * wZ2r
          B_poloidal = BR                * wR2p + BZ                * wZ2p
          vve(i,j,k,n_var + 1) = BR
          vve(i,j,k,n_var + 2) = BZ
          vve(i,j,k,n_var + 3) = Bp
          vve(i,j,k,n_var + 4) = A_radial   
          vve(i,j,k,n_var + 5) = A_poloidal 
          vve(i,j,k,n_var + 6) = U_radial   
          vve(i,j,k,n_var + 7) = U_poloidal 
          vve(i,j,k,n_var + 8) = B_radial   
          vve(i,j,k,n_var + 9) = B_poloidal 
          vve(i,j,k,n_var + 10) = JxB_radial
          vve(i,j,k,n_var + 11) = gradP_radial
          vve(i,j,k,n_var + 12) = RUxB_phi
          vve(i,j,k,n_var + 13) = rho_divU
          vve(i,j,k,n_var + 14) = U_gradRho
#endif            
        end do
        !$omp end parallel do
        
      end do
      
#ifdef fullmhd
#define n_var_max (n_var + n_var_additional)
#else
#define n_var_max (n_var)
#endif
      do l = 1, n_var_max
        call dfftw_plan_dft_r2c_2d(fftw_plan, nequidist_pts, nequidist_tor, vve(:,:,k,l),          &
	  vfour(:,:,k,l), FFTW_ESTIMATE)
! output FFT for nTht/2 +1 poloidal harmonics: check that it is indeed the modes between -nTht/4 and +nTht/4 
        call dfftw_execute(fftw_plan)
        call dfftw_destroy_plan(fftw_plan)
      end do
      
    end do ! k, radial positions

    vfour = vfour / REAL( nequidist_pts * nequidist_tor )
    
    deallocate( vve )
    
  end subroutine transform_qttys
  
  
  
end module fourier

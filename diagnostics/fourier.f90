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
    use mod_model_settings, only: var_A3, var_AR, var_AZ, var_uR, var_uZ, var_up
    use mod_math_operators
    use basis_at_gaussian
    use mod_interp
#endif
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
    real*8 :: F_prof
    real*8 :: ps0_x, ps0_y, psi_abs, wR2r, wZ2r, wR2p, wZ2p
    real*8 :: AR_p, AZ_p, A3_R, A3_Z, AR_Z, AZ_R
    real*8 :: A_radial, A_poloidal, U_radial, U_poloidal, B_radial, B_poloidal
    real*8 :: BR, BZ, Bp
    real*8 :: xjac
#endif

    nequidist_tor = n_plane-1
    nequidist_pts =  nTht
    
    if ( allocated(vfour) ) deallocate(vfour)
#ifdef fullmhd
    call initialise_basis
    allocate(vve(nequidist_pts,nequidist_tor,mapping%nstpts,n_var + 9))
    allocate(vfour(nequidist_pts/2+1,nequidist_tor,mapping%nstpts,n_var + 9))
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
        !$omp     F_prof,                                                                          &
        !$omp     ps0_x, ps0_y, psi_abs, wR2r, wZ2r, wR2p, wZ2p,                                   &
        !$omp     AR_p, AZ_p, A3_R, A3_Z, AR_Z, AZ_R,                                              &
        !$omp     A_radial, A_poloidal, U_radial, U_poloidal, B_radial, B_poloidal,                &
        !$omp     BR, BZ, Bp,                                                                      &
        !$omp     xjac)                                                                            &
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
          call interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
          xjac = op_xjac(R_s, R_t, Z_s, Z_t)
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
#if 1
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
            A3_Z = A3_Z + op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * HZ(iharm,j)
            A3_R = A3_R + op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * HZ(iharm,j)
            call interp(node_list,element_list,i_elm_out,var_AR,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            AR_p = AR_p + eq * HZ_p(iharm,j)
            AR_Z = AR_Z + op_grad_Z(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * HZ(iharm,j)
            call interp(node_list,element_list,i_elm_out,var_AZ,iharm,s_out,t_out,eq,eq_s,eq_t,eq_st,eq_ss,eq_tt)
            AZ_p = AZ_p + eq * HZ_p(iharm,j)
            AZ_R = AZ_R + op_grad_R(eq_s, eq_t, R_s, R_t, Z_s, Z_t, xjac) * HZ(iharm,j)
#endif            
          end do ! iharm
#ifdef fullmhd
          call interp(node_list,element_list,i_elm_out,710,1,s_out,t_out,F_prof,eq_s,eq_t,eq_st,eq_ss,eq_tt)
          BR = ( A3_Z - AZ_p )/ R
          BZ = ( AR_p - A3_R )/ R 
          BP = ( AZ_R - AR_Z )    + F_prof/ R
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
#endif            
        end do
        !$omp end parallel do
        
      end do
      
#ifdef fullmhd
#define n_var_max (n_var + 9)
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

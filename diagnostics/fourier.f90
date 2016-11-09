!> Allows to determine the poloidal magnetic coordinate theta_mag from field line tracing
!! in the n=0 component of the magnetic field and to perform a Fourier transformation of
!! physical quantities in (theta_mag, phi)
module fourier
  
  use constants
  use tr_module 
  use mod_parameters,      only: n_vertex_max, n_order, n_plane, n_tor, n_var, variable_names
  use nodes_elements,  only: node_list, element_list
  use phys_module,     only: F0, xpoint, xcase
  use mod_straight_field_line
  
  implicit none
  
  save
  
  public
  
  integer, parameter   :: FFTW_ESTIMATE = 64  !< (constant of the FFTW library)
  
  
  
  contains
  
  
  
  !> Fourier-transform the physical variables in the magnetic angle theta_mag.
  subroutine transform_qttys(mapping, vfour, m_pol_range)
    
    type(t_theta_mapping), intent(in)    :: mapping
    complex, allocatable,  intent(inout) :: vfour(:,:,:,:) !< Transformed quantities (m,n,irad,ivar)
    integer,               intent(in)    :: m_pol_range(2) !< Range of poloidal mode numbers
    
    real*8, allocatable :: vve(:,:,:,:) ! Variable values (ipol,itor,irad,ivar)
    real*8  :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
    integer :: i, j, k, l, iharm, nn, kv, iv, kf
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail
    real*8  :: v, basis_function
    integer*8 :: fftw_plan
    integer :: nequidist_tor
    
    nequidist_tor = n_plane-1
    nequidist_pts = 2*(maxval(m_pol_range))
    
    if ( allocated(vfour) ) deallocate(vfour)
    allocate(vve(nequidist_pts,nequidist_tor,mapping%nstpts,n_var))
    allocate(vfour(nequidist_pts/2+1,nequidist_tor,mapping%nstpts,n_var))
    vve   = 0.d0
    vfour = 0.d0
    
    
    do k = 1, mapping%nstpts ! radial positions
      
      write(*,'(1x,a,i4)') 'Transforming variables on surface', k
      
      do j = 1, nequidist_tor  ! toroidal positions
        
        !$omp parallel do                                                                          &
	!$omp   default(shared)                                                                    &
	!$omp   firstprivate(i,iharm,R_out,Z_out,i_elm_out,s_out,t_out,ifail,v,l,nn,kv,iv,kf,      &
	!$omp     basis_function,G,G_s,G_t,G_st,G_ss,G_tt)
        do i = 1, nequidist_pts  ! poloidal positions
          
          call find_RZ(node_list,element_list,mapping%rre(k,i-1),mapping%zze(k,i-1),R_out,Z_out,   &
            i_elm_out,s_out,t_out,ifail)
          
          do iharm = 1, n_tor    ! toroidal harmonics

            nn = iharm / 2 ! toroidal mode number (without periodicity)

            if ( iharm == 1 ) then
	      basis_function = 1.d0
            else if ( MOD(iharm,2) == 0 ) then
	      basis_function = cos(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            else
	      basis_function = sin(2.*PI*nn*REAL(j-1)/REAL(nequidist_tor))
            end if
	    
            call basisfunctions2(s_out,t_out,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),   &
              G_ss(1:4,1:4),G_tt(1:4,1:4))
  
	    do l = 1, n_var

              v = 0.d0
  
              do kv = 1, n_vertex_max  ! 4 vertices
    
                iv = element_list%element(i_elm_out)%vertex(kv)  ! the node number
    
                do kf = 1, n_order+1       ! 4 basis functions
    
                  v = v + node_list%node(iv)%values(iharm,kf,l)                                    &
                    * element_list%element(i_elm_out)%size(kv,kf) * G(kv,kf)
    
                end do
    
              end do

              !$omp atomic
              vve(i,j,k,l) = vve(i,j,k,l) + v * basis_function
	    end do
            
          end do
          
        end do
        !$omp end parallel do
        
      end do
      
      do l = 1, n_var
        call dfftw_plan_dft_r2c_2d(fftw_plan, nequidist_pts, nequidist_tor, vve(:,:,k,l),          &
	  vfour(:,:,k,l), FFTW_ESTIMATE)
        call dfftw_execute(fftw_plan)
        call dfftw_destroy_plan(fftw_plan)
      end do
      
    end do

    vfour = vfour / REAL( nequidist_pts * nequidist_tor )
    
    deallocate( vve )
    
  end subroutine transform_qttys
  
  
  
end module fourier

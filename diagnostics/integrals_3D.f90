subroutine Integrals_3D(my_id, node_list,element_list,density_tot,density_in,density_out,pressure,pressure_in,pressure_out)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use constants
use data_structure
use Gauss
use basis_at_gaussian
use phys_module
use pellet_module
use mpi_mod

implicit none


type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8  :: psi_axis, psi_bnd
real*8  :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8  :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8  :: eq_g(n_plane,n_var,n_gauss,n_gauss), eq_s(n_plane,n_var,n_gauss,n_gauss)
real*8  :: eq_t(n_plane,n_var,n_gauss,n_gauss), eq_p(n_plane,n_var,n_gauss,n_gauss)
real*8  :: wgauss_copy(n_gauss)

real*8  :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8  :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8  :: dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz
real*8  :: dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz

integer :: i, j, k, in, ms, mt, mp, iv, inode, ife, n_elements, i_elm_axis, i_elm_xpoint(2), ifail
integer :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
real*8  ::  R_axis,Z_axis,s_axis,t_axis 
real*8  :: current, beta_p, beta_n, beta_t, aminor
real*8  :: xjac, BigR, wst, P_int, C_int, zj0, ps0, r0, T0, T0e, Vol, Volume, Area, Bgeo, psi_limit
real*8  :: density_tot, density_in, density_out,  pressure, pressure_in, pressure_out
real*8  :: current_in, current_out, D_int, D_ext, P_ext, C_ext, P_max, delta_phi, phi, P_tot, D_tot
real*8  :: psi_xpoint(2),R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2)
real*8  :: dTdx, dTdy, drhodx, drhody, dPdx, dPdy, dpsidx, dpsidy
real*8  :: grad_psi, grad_P, grad_P_psi, gradP_psi_max, gradP_max
real*8  :: source_volume, source_pellet
real*8  :: local_pellet_particles, local_plasma_particles, local_pellet_volume  

integer,external :: omp_get_num_threads, omp_get_thread_num


call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI procs

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* Integrals  (3D)                     *'
  write(*,*) '***************************************'
  write(*,*) ' n_plane : ',n_plane
  write(*,*) ' n_cpu   : ',n_cpu
endif

density_tot  = 0.d0
pressure = 0.d0
D_int    = 0.d0
P_int    = 0.d0
C_int    = 0.d0
D_ext    = 0.d0
P_ext    = 0.d0
C_ext    = 0.d0
Vol      = 0.d0
P_tot    = 0.d0
D_tot    = 0.d0
wgauss_copy = wgauss

if (use_pellet) then
  local_pellet_particles = 0.d0
  local_plasma_particles = 0.d0
  local_pellet_volume    = 0.d0
endif

Bgeo = F0 / R_geo

delta_phi = 2.d0 * PI / float(n_plane) / float(n_period)
 
P_max         = 0.d0
gradP_max     = 0.d0
gradP_psi_max = 0.d0

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

if (xpoint) then
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
  psi_bnd  = psi_xpoint(1)
  if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
    psi_bnd = psi_xpoint(2)
  endif
else
  psi_bnd = 0.d0
endif
psi_limit = psi_bnd

ife_delta = ceiling(float(element_list%n_elements) / n_cpu)
ife_min   =      my_id     * ife_delta + 1
ife_max   = min((my_id +1) * ife_delta, element_list%n_elements)

!$omp parallel default(none)                                                                   &
!$omp   shared(element_list,node_list, H, H_s, H_t, HZ, ife_min, ife_max, xpoint, xcase,       &
!$omp          Z_xpoint, my_id, use_pellet, psi_limit, delta_phi, psi_axis, psi_bnd,           &
!$omp          D_tot, D_int, D_Ext, P_tot, P_int, P_ext, Vol, C_int, C_ext,                    &
!$omp          pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi,                       &
!$omp          pellet_radius, pellet_delta_psi, pellet_sig, pellet_length,                     &
!$omp          central_density, pellet_particles,pellet_density, pellet_volume,                &
!$omp          local_pellet_particles, local_plasma_particles, local_pellet_volume,            &
!$omp          wgauss_copy)    &
!$omp   private(ife,iv,inode,element,nodes,i,j, k,in, mp, ms, mt,                              &
!$omp           x_g, y_g, x_s, y_s, x_t, y_t, xjac, eq_g, eq_s, eq_t, eq_p,                    &
!$omp           wst, BigR, r0, T0, T0e, zj0, ps0, dTdx, dTdy, drhodx, drhody, dpsidx, dpsidy,  &
!$omp           dpdx, dpdy, grad_P, grad_psi, grad_P_psi,gradP_max, gradP_psi_max, phi,        &
!$omp           P_max, source_pellet, source_volume, eq_zne, eq_zTe,                           &
!$omp           dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz,    &
!$omp           dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz,    &
!$omp           omp_nthreads,omp_tid)

omp_nthreads = omp_get_num_threads()
omp_tid      = omp_get_thread_num()

!$omp do reduction(+:local_pellet_particles, local_plasma_particles, local_pellet_volume, &
!$omp                D_int, D_ext, P_int, P_ext, C_int, C_ext, Vol, P_tot, D_tot)

do ife = ife_min, ife_max

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1
    
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)

        enddo
      enddo
    enddo
  enddo

  eq_g(:,:,:,:) = 0.d0; eq_s(:,:,:,:) = 0.d0; eq_t(:,:,:,:) = 0.d0; eq_p(:,:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1

      do mp=1,n_plane
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              do in=1,n_tor
                eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
                eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
!                eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
              enddo
            enddo
	    
	  enddo
        enddo
      enddo
      
    enddo
  enddo
  
  do ms=1, n_gauss
    do mt=1, n_gauss
 
      call density(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                   dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

      call temperature(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
  
    enddo
  enddo

  eq_zTe = eq_zTe / 2.d0	! electron temperature	

!--------------------------------------------------- sum over the Gaussian integration points
  
  do mp=1,n_plane
    
    phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)

    do ms=1, n_gauss

      do mt=1, n_gauss

        wst = wgauss_copy(ms)*wgauss_copy(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        r0     = eq_g(mp,5,ms,mt)
        T0     = eq_g(mp,6,ms,mt) 
        T0e    = eq_g(mp,6,ms,mt) /2.d0
        zj0    = eq_g(mp,3,ms,mt)
        ps0    = eq_g(mp,1,ms,mt)
	
        dTdx   = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
        dTdy   = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
        drhodx = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
        drhody = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac

        dpsidx = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
        dpsidy = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
	
	dPdx = r0 * dTdx + T0 * drhodx
	dPdy = r0 * dTdy + T0 * drhody
	
	grad_P   = sqrt(dPdx**2   + dPdy**2) 
	grad_psi = sqrt(dpsidx**2 + dpsidy**2)
	
	grad_P_psi = (dPdx * dpsidx + dPdy * dpsidy)/grad_psi

        P_tot = P_tot + r0 * T0 * xjac * BigR * wst * delta_phi
        D_tot = D_tot + r0      * xjac * BigR * wst * delta_phi
      
        P_max = max(P_max,r0 * T0)
	
	gradP_max     = max(gradP_max,grad_P)
	gradP_psi_max = max(gradP_psi_max,grad_P_psi)    
        
        if (use_pellet) then
        
!          call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
!                              pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
!                              x_g(ms,mt),y_g(ms,mt), ps0, phi, r0, T0e, central_density, &
!                              pellet_particles, pellet_density, pellet_volume, source_pellet, source_volume)                  

          call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                              pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
                              x_g(ms,mt),y_g(ms,mt), ps0, phi, eq_zne(ms,mt), eq_zTe(ms,mt), central_density, &
                              pellet_particles, pellet_density, pellet_volume, source_pellet, source_volume)                  

          local_pellet_particles = local_pellet_particles + source_pellet * bigR * xjac * wst * delta_phi
          local_plasma_particles = local_plasma_particles + r0            * bigR * xjac * wst * delta_phi
          local_pellet_volume    = local_pellet_volume    + source_volume * bigR * xjac * wst * delta_phi 

        endif

        if (ps0 .lt. psi_limit) then
        
          D_int = D_int + r0        * xjac * BigR * wst * delta_phi
          P_int = P_int + r0 * T0   * xjac * BigR * wst * delta_phi
          C_int = C_int + zj0 /BigR * xjac * BigR * wst * delta_phi                
          Vol   = Vol   +             xjac * BigR * wst * delta_phi
        
        else

          D_ext = D_ext + r0         * xjac * BigR * wst * delta_phi     
          P_ext = P_ext + r0   * T0  * xjac * BigR * wst * delta_phi
          C_ext = C_ext + zj0 / BigR * xjac * BigR * wst * delta_phi
        
        endif
      
      enddo
    enddo
  enddo
    
enddo
!$omp end do
!$omp end parallel

call MPI_AllReduce(D_int,density_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(D_ext,density_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_int,pressure_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_ext,pressure_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(C_int,current_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(C_ext,current_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(Vol,Volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(D_tot,density_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_tot,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

if (use_pellet) then
  call MPI_AllReduce(local_pellet_particles,total_pellet_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_plasma_particles,total_plasma_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_pellet_volume,total_pellet_volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
endif

current = C_int / MU_zero

if (my_id .eq. 0) then
  write(*,'(A,2e14.6)') ' Integrals_3D, PELLET : ',pellet_volume, total_pellet_volume
  if (index_start .gt.0) then
    write(*,'(A,8e14.6)') ' Volume   : ',xtime(index_start),volume
    write(*,'(A,8e14.6)') ' density  (total/in/out) : ',xtime(index_start),density_tot,  density_in,  density_out 
    write(*,'(A,8e14.6)') ' pressure (total/in/out) : ',xtime(index_start),pressure, pressure_in, pressure_out, P_max, gradP_max, gradP_psi_max
    write(*,'(A,8e14.6)') ' current  (in/out)       : ',xtime(index_start),current_in, current_out 
 else 
    write(*,'(A,8e14.6)') ' Volume   : ',0.d0,volume
    write(*,'(A,8e14.6)') ' density  (total/in/out) : ',0.d0,density_tot,  density_in,  density_out 
    write(*,'(A,8e14.6)') ' pressure (total/in/out) : ',0.d0,pressure, pressure_in, pressure_out, P_max, gradP_max, gradP_psi_max
    write(*,'(A,8e14.6)') ' current  (in/out)       : ',0.d0,current_in, current_out 
  endif
endif
return
end

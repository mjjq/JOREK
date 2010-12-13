subroutine Integrals_3D(node_list,element_list,&
                     density,density_in,density_out,pressure,pressure_in,pressure_out)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use Gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: eq_g(n_plane,n_var,n_gauss,n_gauss), eq_s(n_plane,n_var,n_gauss,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss,n_gauss), eq_p(n_plane,n_var,n_gauss,n_gauss)

integer    :: i, j, k, in, ms, mt, mp, iv, inode, ife, n_elements, i_elm_xpoint, ifail
real*8     :: current, beta_p, beta_n, beta_t, MU_zero, aminor
real*8     :: xjac, BigR, wst, P_int, C_int, ZJ_0, PS_0, Volume, Area, PI, Bgeo, psi_limit
real*8     :: density, density_in, density_out,  pressure, pressure_in, pressure_out
real*8     :: current_in, current_out, D_int, D_ext, P_ext, C_ext, P_max, delta_phi
real*8     :: psi_xpoint,R_xpoint,Z_xpoint,s_xpoint,t_xpoint
real*8     :: dTdx, dTdy, drhodx, drhody, dPdx, dPdy, dpsidx, dpsidy
real*8     :: grad_psi, grad_P, grad_P_psi, gradP_psi_max, gradP_max

write(*,*) '***************************************'
write(*,*) '* Integrals  (3D)                     *'
write(*,*) '***************************************'
write(*,*) ' n_plane : ',n_plane


density  = 0.d0
pressure = 0.d0
D_int    = 0.d0
P_int    = 0.d0
C_int    = 0.d0
D_ext    = 0.d0
P_ext    = 0.d0
C_ext    = 0.d0
Volume   = 0.d0

Bgeo = F0 / R_geo

PI = 2.d0*asin(1.d0)
 
delta_phi = 2.d0 * PI / float(n_plane)
 
P_max     = 0.d0
gradP_max = 0.d0
gradP_psi_max = 0.d0

if (xpoint) then
  call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
  psi_limit = psi_xpoint
else
  psi_limit = 0.d0
endif
write(*,*) ' psi_limit : ',psi_limit

do ife =1,  element_list%n_elements

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
!--------------------------------------------------- sum over the Gaussian integration points
  
  do mp=1,n_plane

    do ms=1, n_gauss

      do mt=1, n_gauss

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        rho_0 = eq_g(mp,5,ms,mt)
        T_0   = eq_g(mp,6,ms,mt)
        ZJ_0  = eq_g(mp,3,ms,mt)
        PS_0  = eq_g(mp,1,ms,mt)
	
        dTdx   = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
        dTdy   = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
        drhodx = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
        drhody = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac

        dpsidx = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
        dpsidy = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
	
	dPdx = rho_0 * dTdx + T_0 * drhodx
	dPdy = rho_0 * dTdy + T_0 * drhody
	
	grad_P   = sqrt(dPdx**2   + dPdy**2) 
	grad_psi = sqrt(dpsidx**2 + dpsidy**2)
	
	grad_P_psi = (dPdx * dpsidx + dPdy * dpsidy)/grad_psi

        pressure = pressure + rho_0 * T_0 * xjac * BigR * wst * delta_phi
        density  = density  + rho_0       * xjac * BigR * wst * delta_phi
      
        P_max = max(P_max,rho_0 * T_0)
	
	gradP_max     = max(gradP_max,grad_P)
	gradP_psi_max = max(gradP_psi_max,grad_P_psi)

        if (PS_0 .lt. psi_limit) then

          D_int = D_int + rho_0       * xjac * BigR * wst * delta_phi
          P_int = P_int + rho_0 * T_0 * xjac * BigR * wst * delta_phi
          C_int = C_int + ZJ_0 /BigR  * xjac * BigR * wst * delta_phi
                
          Volume = Volume + BigR * xjac * wst * delta_phi
        
        else

          D_ext = D_ext + rho_0       * xjac * BigR * wst * delta_phi     
          P_ext = P_ext + rho_0 * T_0 * xjac * BigR * wst * delta_phi
          C_ext = C_ext + ZJ_0 /BigR  * xjac * BigR * wst * delta_phi
        
        endif
      
      enddo
    enddo
  enddo
    
enddo

density_in   = D_int
density_out  = D_ext
pressure_in  = P_int
pressure_out = P_ext
current_in   = C_int
current_out  = C_ext

MU_zero = 4.d0*PI * 1.d-7

current = C_int / MU_zero

if (index_start .gt.0) then
  write(*,'(A,8e14.6)') ' Volume   : ',xtime(index_start),volume
  write(*,'(A,8e14.6)') ' density  (total/in/out) : ',xtime(index_start),density,  density_in,  density_out 
  write(*,'(A,8e14.6)') ' pressure (total/in/out) : ',xtime(index_start),pressure, pressure_in, pressure_out, P_max, gradP_max, gradP_psi_max
  write(*,'(A,8e14.6)') ' current  (in/out)       : ',xtime(index_start),current_in, current_out 
else
  write(*,'(A,8e14.6)') ' Volume   : ',0.d0,volume
  write(*,'(A,8e14.6)') ' density  (total/in/out) : ',0.d0,density,  density_in,  density_out 
  write(*,'(A,8e14.6)') ' pressure (total/in/out) : ',0.d0,pressure, pressure_in, pressure_out, P_max, gradP_max, gradP_psi_max
  write(*,'(A,8e14.6)') ' current  (in/out)       : ',0.d0,current_in, current_out 
endif
return
end

!> Determines some integrals over the JOREK computational domain to determine the total current etc.
subroutine integrals(node_list, element_list, psi_axis, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  aminor, Bgeo, current, beta_p, beta_t, beta_n, density, density_in, density_out, pressure,       &
  pressure_in, pressure_out)

use constants
use parameters
use data_structure
use Gauss
use basis_at_gaussian
use phys_module
use domains

implicit none

! --- Routine parameters
type(type_node_list),    intent(in)    :: node_list
type(type_element_list), intent(in)    :: element_list
real*8,                  intent(in)    :: psi_axis
real*8,                  intent(in)    :: R_xpoint(2)
real*8,                  intent(in)    :: Z_xpoint(2)
real*8,                  intent(in)    :: psi_xpoint(2)
real*8,                  intent(in)    :: psi_limit
real*8,                  intent(in)    :: aminor
real*8,                  intent(out)   :: Bgeo
real*8,                  intent(out)   :: current
real*8,                  intent(out)   :: beta_p
real*8,                  intent(out)   :: beta_t
real*8,                  intent(out)   :: beta_n
real*8,                  intent(out)   :: density
real*8,                  intent(out)   :: density_in
real*8,                  intent(out)   :: density_out
real*8,                  intent(out)   :: pressure
real*8,                  intent(out)   :: pressure_in
real*8,                  intent(out)   :: pressure_out

! --- Local variables
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: eq_g(n_var,n_gauss,n_gauss), eq_s(n_var,n_gauss,n_gauss), eq_t(n_var,n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: xjac, BigR, wst, P_int, C_intern, ZJ_0, PS_0, Volume, Area
real*8     :: rho_00, T_00, Ti_00, Te_00, current_in, current_out 
real*8     :: C_hel, P_hel, D_int, D_ext, P_ext, C_ext
real*8     :: particle_source, heat_source, heating_power

write(*,*) '***************************************'
write(*,*) '* Integrals                           *'
write(*,*) '***************************************'

density  = 0.d0
pressure = 0.d0
D_int    = 0.d0
P_int    = 0.d0
C_intern = 0.d0
D_ext    = 0.d0
P_ext    = 0.d0
C_ext    = 0.d0
P_hel    = 0.d0
C_hel    = 0.d0
Volume   = 0.d0
Area     = 0.d0
heating_power = 0

Bgeo = F0 / R_geo

do ife =1, element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

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

  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss

          do k=1,n_var
            eq_g(k,ms,mt)  = eq_g(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H(i,j,ms,mt)
            eq_s(k,ms,mt)  = eq_s(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
            eq_t(k,ms,mt)  = eq_t(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H_t(i,j,ms,mt)
          enddo

        enddo
      enddo
    enddo
  enddo
!--------------------------------------------------- sum over the Gaussian integration points

  do ms=1, n_gauss

    do mt=1, n_gauss

      wst = wgauss(ms)*wgauss(mt)

      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      BigR = x_g(ms,mt)

      rho_00 = eq_g(5,ms,mt)
      if (jorek_model .eq. 400) then
        Ti_00 = eq_g(6,ms,mt)
        Te_00 = eq_g(n_var,ms,mt)
        T_00  = Ti_00 + Te_00
      else
        T_00  = eq_g(6,ms,mt)
      endif
      ZJ_0  = eq_g(3,ms,mt)
      PS_0  = eq_g(1,ms,mt)
      
      pressure = pressure + rho_00 * T_00 * xjac * BigR * wst
      density  = density  + rho_00       * xjac * BigR * wst
      
      if ( in_plasma(x_g(ms,mt),y_g(ms,mt),eq_g(1,ms,mt),xpoint,&
        xcase,R_xpoint,Z_xpoint,psi_xpoint,psi_limit) ) then
        
        call sources(xpoint, xcase, eq_g(2,ms,mt), Z_xpoint, eq_g(1,ms,mt), psi_axis, &
          psi_limit, particle_source, heat_source)
        
        ! --- 3D integrals
        D_int = D_int + rho_00       * xjac * BigR * wst
        P_int = P_int + rho_00 * T_00 * xjac * BigR * wst
        C_intern = C_intern + ZJ_0 /BigR  * xjac * BigR * wst
        
        ! --- 2D integrals
        P_hel = P_hel + rho_00 * T_00 * xjac * wst
        C_hel = C_hel + ZJ_0 /BigR  * xjac * wst
        Volume = Volume + 2.d0 * PI * BigR * xjac * wst
        heating_power = heating_power + 2.d0 * PI * BigR * xjac * wst * heat_source
        Area   = Area   + xjac * wst
        
      else
        
        D_ext = D_ext + rho_00       * xjac * BigR * wst      
        P_ext = P_ext + rho_00 * T_00 * xjac * BigR * wst
        C_ext = C_ext + ZJ_0 /BigR  * xjac * BigR * wst
        
      endif
      
    enddo
  enddo
enddo

density      = density  * 2.d0 * PI
density_in   = D_int    * 2.d0 * PI
density_out  = D_ext    * 2.d0 * PI
pressure     = pressure * 2.d0 * PI
pressure_in  = P_int    * 2.d0 * PI
pressure_out = P_ext    * 2.d0 * PI
current_in   = C_intern * 2.d0 * PI
current_out  = C_ext    * 2.d0 * PI

current = C_hel / MU_zero
beta_p  = 8.d0 * PI * P_hel / (C_hel**2 )
beta_t  = 2.d0 * P_hel / Bgeo**2 / (Area)
beta_n  = 100.d0 * (4.*PI/10.) * beta_t / (MU_zero * abs(current) /  (aminor * Bgeo))

write(*,'(A,f12.7)') ' psi_limit: ',psi_limit
write(*,'(A,f12.7,A)') ' current  : ',current/1.e6,' MA'
write(*,'(A,f12.7)') ' beta_p   : ',beta_p
write(*,'(A,f12.7)') ' beta_t   : ',beta_t
write(*,'(A,f12.7,A)') ' beta_n   : ',beta_n,' [%]'
write(*,'(A,f12.7,A)') ' Area     : ',area,' m^2'
write(*,'(A,f12.7,A)') ' Volume   : ',volume,' m^3'
write(*,'(A,es18.7,A)') ' Heating power : ',heating_power,' / sqrt((mu_0)^3 rho_0) W'

write(*,'(A,5f10.5)') ' density  (total/in/out) : ',density,  density_in,  density_out 
write(*,'(A,5f10.5)') ' pressure (total/in/out) : ',pressure, pressure_in, pressure_out 
write(*,'(A,5f10.5)') ' current  (in/out)       : ',current_in, current_out 

return
end subroutine integrals

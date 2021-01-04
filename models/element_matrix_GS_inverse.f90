subroutine element_matrix_GS_inverse(xpoint,xcase,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module, only: treat_axis

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: factor(n_gauss,n_gauss)
real*8     :: eq_g(n_gauss,n_gauss),   eq_s(n_gauss,n_gauss),   eq_t(n_gauss,n_gauss)
real*8     :: eq2_g(n_gauss,n_gauss),  eq2_s(n_gauss,n_gauss),  eq2_t(n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1))

real*8     :: xjac, wst
real*8     :: ps0_x, ps0_y, v, v_x, v_y, psi, psi_x, psi_y, rhs_ij
integer    :: ms, mt, i, j, k, l, index_ij, index_kl, itype, ivar_in, ivar_out, i_harm, xcase
logical    :: xpoint
real*8     :: Z_xpoint(2),psi_axis,psi_bnd,dj_dpsi,dj_dz
real*8     :: zn,dn_dpsi,dn_dz, zT,dT_dpsi,dT_dz, zFFprime,dFFprime_dpsi,dFFprime_dz

! axis treatment related variables
integer    :: inode
real*8     :: esize(n_vertex_max,n_order+1)
real*8     :: BasFun   (n_vertex_max, n_order+1, n_gauss, n_gauss)
real*8     :: BasFun_s (n_vertex_max, n_order+1, n_gauss, n_gauss)
real*8     :: BasFun_t (n_vertex_max, n_order+1, n_gauss, n_gauss)

ELM=0.d0
RHS=0.d0

! change basis function for elements on the grid axis
esize(:,:) = element%size(:,:)
BasFun  = H ; BasFun_s  = H_s ; BasFun_t  = H_t

if(treat_axis .and. element%axis_element)then
  call on_the_axis(element, nodes, H  ,  BasFun  )
  call on_the_axis(element, nodes, H_s,  BasFun_s)
  call on_the_axis(element, nodes, H_t,  BasFun_t)
  esize(1  ,:) = 1.0d0
  esize(2:3,:) = element%size(2:3,:)
  esize(4  ,:) = 1.0d0
endif

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g(:,:)   = 0.d0; x_s(:,:)   = 0.d0; x_t(:,:)   = 0.d0;
y_g(:,:)   = 0.d0; y_s(:,:)   = 0.d0; y_t(:,:)   = 0.d0;
eq_g(:,:)  = 0.d0; eq_s(:,:)  = 0.d0; eq_t(:,:)  = 0.d0;
eq2_g(:,:) = 0.d0; eq2_s(:,:) = 0.d0; eq2_t(:,:) = 0.d0;

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

       ! change basis function for axis elements only for physical variables
       eq_g(ms,mt)  = eq_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * esize(i,j) * BasFun(i,j,ms,mt)
       eq_s(ms,mt)  = eq_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * esize(i,j) * BasFun_s(i,j,ms,mt)
       eq_t(ms,mt)  = eq_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * esize(i,j) * BasFun_t(i,j,ms,mt)

       eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * esize(i,j) * BasFun(i,j,ms,mt)
       eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * esize(i,j) * BasFun_s(i,j,ms,mt)
       eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * esize(i,j) * BasFun_t(i,j,ms,mt)

     enddo
   enddo
 enddo
enddo


factor = 1.d0/ x_g                    ! Grad-Shafranov

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

   ps0_x = (   y_t(ms,mt) * eq_s(ms,mt) - y_s(ms,mt) * eq_t(ms,mt) ) / xjac
   ps0_y = ( - x_t(ms,mt) * eq_s(ms,mt) + x_s(ms,mt) * eq_t(ms,mt) ) / xjac

   do i=1,n_vertex_max

     do j=1,n_order+1

       index_ij = (i-1)*(n_order+1) + j

       ! change basis function for axis elements
       v   = BasFun(i,j,ms,mt)  * esize(i,j)
       v_x = (  y_t(ms,mt) * BasFun_s(i,j,ms,mt) - y_s(ms,mt) * BasFun_t(i,j,ms,mt) ) * esize(i,j) / xjac
       v_y = (- x_t(ms,mt) * BasFun_s(i,j,ms,mt) + x_s(ms,mt) * BasFun_t(i,j,ms,mt) ) * esize(i,j) / xjac

       RHS(index_ij) = RHS(index_ij) + (v_x * ps0_x + v_y * ps0_y) * factor(ms,mt) * xjac * wst
       RHS(index_ij) = RHS(index_ij) + v * eq2_g(ms,mt) * factor(ms,mt) * xjac * wst    ! solve for perturbation only

       do k=1,n_vertex_max

         do l=1,n_order+1

           ! change basis function for axis elements
           psi   = BasFun(k,l,ms,mt)  * esize(k,l)
           psi_x = (   y_t(ms,mt) * BasFun_s(k,l,ms,mt) - y_s(ms,mt) * BasFun_t(k,l,ms,mt) ) * esize(k,l) / xjac
           psi_y = ( - x_t(ms,mt) * BasFun_s(k,l,ms,mt) + x_s(ms,mt) * BasFun_t(k,l,ms,mt) ) * esize(k,l) / xjac

           index_kl = (k-1)*(n_order+1) + l

           ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) - psi * v * factor(ms,mt) * xjac * wst

         enddo
       enddo

     enddo
   enddo

 enddo
enddo

return
end

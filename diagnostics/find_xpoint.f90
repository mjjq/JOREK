subroutine find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: psi_xpoint, dpsi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint, ps_s, ps_t
real*8  :: grad_psi, grad_psi_min
real*8  :: R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: ps_x, ps_y, xjac
integer :: i_elm_xpoint, ij_xpoint(2), i, iv, ms, mt, kf, kv, ifail

real*8  :: x(2), s, t, xerr, ferr, rs_tolerance

logical :: early_exit
parameter (rs_tolerance = 1.d-8)

write(*,*) '*********************************'
write(*,*) '*     find_xpoint               *'
write(*,*) '*********************************'


dpsi_xpoint  = 1.d20
grad_psi_min = 1.d20
Z_xpoint     = 0.d0

do i=1,element_list%n_elements

  do ms = 1, 4           ! 4 Gaussian points
    do mt = 1, 4         ! 4 Gaussian points

      ps_s = 0.d0
      ps_t = 0.d0
      R_s  = 0.d0 
      Z_s  = 0.d0
      R_t  = 0.d0 
      Z_t  = 0.d0
      Z    = 0.d0

      do kf = 1, 4       ! 4 basis functions
        do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          ps_s = ps_s + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          ps_t = ps_t + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

          Z   = Z   + node_list%node(iv)%x(kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

          R_s = R_s + node_list%node(iv)%x(kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          Z_s = Z_s + node_list%node(iv)%x(kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          R_t = R_t + node_list%node(iv)%x(kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
          Z_t = Z_t + node_list%node(iv)%x(kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

        enddo
      enddo

      xjac = R_s * Z_t - R_t * Z_s
      ps_x = (  ps_s * Z_t - ps_t * Z_s)/ xjac
      ps_y = (- ps_s * R_t + ps_t * R_s)/ xjac

      grad_psi = sqrt(ps_x*ps_x + ps_y*ps_y)
      
      if ((grad_psi .lt. grad_psi_min) .and. ((Z) .lt. -0.4d0)) then
        grad_psi_min = grad_psi
	Z_xpoint     = Z
        i_elm_xpoint = i
        ij_xpoint(1) = ms;         ij_xpoint(2)  = mt
      endif

    enddo
  enddo
  
enddo

s=Xgauss(ij_xpoint(1)) ; t=Xgauss(ij_xpoint(2))

call mnewtax(node_list,element_list,i_elm_xpoint,s,t,xerr,ferr,ifail)

if (ifail .ne. 0 ) write(*,*) ' MNEWTAX : ifail = ',ifail

call interp(node_list,element_list,i_elm_xpoint,1,1,s,t,psi_xpoint,P_s,P_t,P_st,P_ss,P_tt)

call interp_RZ(node_list,element_list,i_elm_xpoint,s,t,R_xpoint,R_s,R_t,R_st,R_ss,R_tt,Z_xpoint,Z_s,Z_t,Z_st,Z_ss,Z_tt)

s_xpoint = s
t_xpoint = t

xjac = R_s * Z_t - R_t * Z_s
ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
ps_y = (- P_s * R_t + P_t * R_s)/ xjac

write(*,'(A,i6,4f14.8)') ' X-point : ',i_elm_xpoint,R_xpoint,Z_xpoint,psi_xpoint,sqrt(ps_x**2+ps_y**2)

if (sqrt(ps_x**2+ps_y**2) .gt. 1.d-4) ifail=1

if (ifail .ne. 0 ) write(*,*) ' find_xpoint : ifail = ',ifail

return
END

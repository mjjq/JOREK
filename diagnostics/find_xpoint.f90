!> Routine determines the position(s) of the xpoint(s).
subroutine find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)

use data_structure
use gauss
use basis_at_gaussian

implicit none

! --- Routine parameters
integer,                  intent(in)    :: my_id
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(out)   :: psi_xpoint(2)
real*8,                   intent(out)   :: R_xpoint(2)
real*8,                   intent(out)   :: Z_xpoint(2)
integer,                  intent(out)   :: i_elm_xpoint(2)
real*8,                   intent(out)   :: s_xpoint(2)
real*8,                   intent(out)   :: t_xpoint(2)
integer,                  intent(in)    :: xcase
integer,                  intent(out)   :: ifail

! --- Local variables
real*8  :: grad_psi, grad_psi_min(2), ps_s, ps_t
real*8  :: R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: ps_x, ps_y, xjac
integer :: ij_xpoint(2,2), i, iv, ms, mt, kf, kv

real*8  :: x(2), s, t, xerr, ferr
logical :: early_exit

real*8, parameter :: rs_tolerance = 1.d-8

if (my_id .eq. 0) then
  write(*,*) '*********************************'
  write(*,*) '*     find_xpoint               *'
  write(*,*) '*********************************'
endif

ifail = 0

grad_psi_min = 1.d20
Z_xpoint(1)  = 0.d0
Z_xpoint(2)  = 0.d0

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
      
      ! --- Look for the lower Xpoint
      if ((grad_psi .lt. grad_psi_min(1)) .and. (Z .lt. -0.4d0) .and. (xcase .ne. 2)) then
        grad_psi_min(1) = grad_psi
	Z_xpoint(1)     = Z
        i_elm_xpoint(1) = i
        ij_xpoint(1,1) = ms;         ij_xpoint(1,2)  = mt
      endif
      ! --- And for the upper Xpoint
      if ((grad_psi .lt. grad_psi_min(2)) .and. (Z .gt.  0.4d0) .and. (xcase .ne. 1)) then
        grad_psi_min(2) = grad_psi
	Z_xpoint(2)     = Z
        i_elm_xpoint(2) = i
        ij_xpoint(2,1) = ms;         ij_xpoint(2,2)  = mt
      endif

    enddo
  enddo
  
enddo

if(xcase .ne. 2) then
  s=Xgauss(ij_xpoint(1,1)) ; t=Xgauss(ij_xpoint(1,2))
  call mnewtax(node_list,element_list,i_elm_xpoint(1),s,t,xerr,ferr,ifail)
  if ((ifail .ne. 0 ).and.(my_id .eq.0)) write(*,*) ' MNEWTAX LowerXpoint: ifail = ',ifail

  call interp(node_list,element_list,i_elm_xpoint(1),1,1,s,t,psi_xpoint(1),P_s,P_t,P_st,P_ss,P_tt)
  call interp_RZ(node_list,element_list,i_elm_xpoint(1),s,t,R_xpoint(1),R_s,R_t,R_st,R_ss,R_tt,Z_xpoint(1),Z_s,Z_t,Z_st,Z_ss,Z_tt)
  s_xpoint(1) = s
  t_xpoint(1) = t
  
  xjac = R_s * Z_t - R_t * Z_s
  ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
  ps_y = (- P_s * R_t + P_t * R_s)/ xjac
  
  if (my_id .eq. 0) then
    write(*,'(A,i6,4f14.8)') ' Lower X-point : ',i_elm_xpoint(1),R_xpoint(1),Z_xpoint(1),psi_xpoint(1),sqrt(ps_x**2+ps_y**2)
  endif
  if (sqrt(ps_x**2+ps_y**2) .gt. 1.d-4) ifail=1
  if ((ifail .ne. 0 ).and.(my_id .eq.0)) write(*,*) ' find_xpoint : LowerXpoint ifail = ',ifail
endif

if(xcase .ne. 1) then
  s=Xgauss(ij_xpoint(2,1)) ; t=Xgauss(ij_xpoint(2,2))
  call mnewtax(node_list,element_list,i_elm_xpoint(2),s,t,xerr,ferr,ifail)
  if ((ifail .ne. 0 ).and.(my_id .eq.0)) write(*,*) ' MNEWTAX UpperXpoint: ifail = ',ifail

  call interp(node_list,element_list,i_elm_xpoint(2),1,1,s,t,psi_xpoint(2),P_s,P_t,P_st,P_ss,P_tt)
  call interp_RZ(node_list,element_list,i_elm_xpoint(2),s,t,R_xpoint(2),R_s,R_t,R_st,R_ss,R_tt,Z_xpoint(2),Z_s,Z_t,Z_st,Z_ss,Z_tt)
  s_xpoint(2) = s
  t_xpoint(2) = t
  
  xjac = R_s * Z_t - R_t * Z_s
  ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
  ps_y = (- P_s * R_t + P_t * R_s)/ xjac
  
  if (my_id .eq. 0) then
    write(*,'(A,i6,4f14.8)') ' Upper X-point : ',i_elm_xpoint(2),R_xpoint(2),Z_xpoint(2),psi_xpoint(2),sqrt(ps_x**2+ps_y**2)
  endif
  if (sqrt(ps_x**2+ps_y**2) .gt. 1.d-4) ifail=1
  if ((ifail .ne. 0 ).and.(my_id .eq.0)) write(*,*) ' find_xpoint : UpperXpoint ifail = ',ifail
endif


return
end subroutine find_xpoint

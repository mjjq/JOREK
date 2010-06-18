subroutine find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: psi_xpoint, dpsi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint, ps_s, ps_t
real*8  :: DPS(n_order+1,n_order+1), dpsi_min,  Z(n_order+1,n_order+1), Z_min
real*8  :: R_s, R_t, R_st, R_ss, R_tt, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
integer :: i_elm_xpoint, ij_xpoint(2), i, iv, ms, mt, kf, kv, ifail

real*8  :: x(2), s, t, xerr, ferr, rs_tolerance

logical :: early_exit
parameter (rs_tolerance = 1.d-8)

write(*,*) '*********************************'
write(*,*) '*     find_xpoint               *'
write(*,*) '*********************************'


dpsi_xpoint = 1.d20
                                                                                                    
ij_xpoint = 1 ! XL: In some cases it is not initialised... 
i_elm_xpoint = 1 ! XL : In some cases it is not initialised... 
do i=1,element_list%n_elements

  Z = 0.d0

  do ms = 1, 4           ! 4 Gaussian points
    do mt = 1, 4         ! 4 Gaussian points

      ps_s = 0.d0
      ps_t = 0.d0

      do kf = 1, 4       ! 4 basis functions
        do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          ps_s = ps_s + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          ps_t = ps_t + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

          Z(ms,mt) = Z(ms,mt)  + node_list%node(iv)%x(kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

        enddo
      enddo

      DPS(ms,mt) = ps_s*ps_s + ps_t*ps_t

    enddo
  enddo

  dpsi_min = minval(DPS)
  Z_min    = minval(Z)

  if ((dpsi_min .lt. dpsi_xpoint) .and. (Z_min .lt. -0.8d0))  then
    dpsi_xpoint   = dpsi_min
    i_elm_xpoint  = i
    ij_xpoint     = minloc(DPS)
    Z_xpoint      = Z_min
  endif

enddo

!write(*,*) ' estimate : ',Z_xpoint

s=Xgauss(ij_xpoint(1)) ; t=Xgauss(ij_xpoint(2))
call mnewtax(node_list,element_list,i_elm_xpoint,s,t,xerr,ferr,ifail)

if (ifail .ne. 0 ) write(*,*) ' MNEWTAX : ifail = ',ifail

call interp(node_list,element_list,i_elm_xpoint,1,1,s,t,psi_xpoint,P_s,P_t,P_st,P_ss,P_tt)

call interp_RZ(node_list,element_list,i_elm_xpoint,s,t,R_xpoint,R_s,R_t,R_st,R_ss,R_tt,Z_xpoint,Z_s,Z_t,Z_st,Z_ss,Z_tt)

s_xpoint = s
t_xpoint = t

write(*,'(A,i6,4f14.8)') ' X-point : ',i_elm_xpoint,R_xpoint,Z_xpoint,psi_xpoint

return
END

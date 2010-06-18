subroutine find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian

implicit none

interface
   subroutine mnewtax(node_list,element_list,i_elm, r, s, errx, errf, ifail)
     !-------------------------------------------------------------------------
     ! solves two non-linear equations using Newtons method (from numerical recipes)
     ! LU decomposition replaced by explicit solution of 2x2 matrix.
     !-------------------------------------------------------------------------
     use data_structure
     
     type (type_node_list)    :: node_list
     type (type_element_list) :: element_list

     real*8    :: r, s
     real*8    :: errf, errx
     integer   :: ifail, i_elm
   end subroutine mnewtax
end interface

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
real*8  :: PS(n_order+1,n_order+1), psi_min, psi_max, Y(n_order+1,n_order+1), Y_min
real*8  :: R_s, R_t, R_st, R_ss, R_tt, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
integer :: i_elm_axis, ij_axis(2), i, iv, ms, mt, kf, kv, ifail

real*8  :: x(2), s, t, xerr, ferr, rs_tolerance

logical :: early_exit
parameter (rs_tolerance = 1.d-8)

write(*,*) '*********************************'
write(*,*) '*     find_axis                 *'
write(*,*) '*********************************'

ij_axis = 1 ! XL: In some cases it is not initialised...
psi_axis = 1.d20

do i=1,element_list%n_elements

  PS = 0.d0
  Y  = 0.d0

  do ms = 1, 4           ! 4 Gaussian points
    do mt = 1, 4         ! 4 Gaussian points

      do kf = 1, 4       ! 4 basis functions
        do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          PS(ms,mt)  = PS(ms,mt) + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          Y(ms,mt)   = Y(ms,mt)  + node_list%node(iv)%x(kf,2)         * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

        enddo
      enddo

    enddo
  enddo

  psi_min = minval(PS)
  Y_min   = minval(Y)

  if ((psi_min .lt. psi_axis) .and. (y_min .ge. -0.5))  then
    psi_axis    = psi_min
    i_elm_axis  = i
    ij_axis     = minloc(PS)
  endif
!write(*,'(A,4f14.8)') ' magnetic axisdd : ',R_axis,Z_axis,psi_axis
enddo

s=Xgauss(ij_axis(1)) ; t=Xgauss(ij_axis(2))

call mnewtax(node_list,element_list,i_elm_axis,s,t,xerr,ferr,ifail)

if (ifail .ne. 0 ) write(*,*) ' MNEWTAX : ifail = ',ifail

call interp(node_list,element_list,i_elm_axis,1,1,s,t,psi_axis,P_s,P_t,P_st,P_ss,P_tt)

call interp_RZ(node_list,element_list,i_elm_axis,s,t,R_axis,R_s,R_t,R_st,R_ss,R_tt,Z_axis,Z_s,Z_t,Z_st,Z_ss,Z_tt)

s_axis = s
t_axis = t

write(*,'(A,4f14.8)') ' magnetic axis : ',R_axis,Z_axis,psi_axis

return
END

!> Locate the position of the magnetic axis.
subroutine find_axis(my_id, node_list, element_list, psi_axis, R_axis, Z_axis, i_elm_axis, s_axis, &
  t_axis, ifail)

use data_structure
use gauss
use basis_at_gaussian
use phys_module, only: R_geo, tokamak_device

implicit none

interface
   subroutine mnewtax(node_list,element_list,i_elm, r, s, errx, errf, ifail)
     use data_structure
     type (type_node_list)    :: node_list
     type (type_element_list) :: element_list
     real*8    :: r, s, errf, errx
     integer   :: ifail, i_elm
   end subroutine mnewtax
end interface

! --- Routine parameters
integer,                 intent(in)  :: my_id        !< MPI proc number
type(type_node_list),    intent(in)  :: node_list    !< List of grid nodes
type(type_element_list), intent(in)  :: element_list !< List of grid elements
real*8,                  intent(out) :: psi_axis     !< Poloidal flux at axis
real*8,                  intent(out) :: R_axis       !< R-position of axis
real*8,                  intent(out) :: Z_axis       !< Z-position of axis
real*8,                  intent(out) :: s_axis       !< s-position of axis in Bezier element i_elm_axis
real*8,                  intent(out) :: t_axis       !< t-position of axis in Bezier element i_elm_axis
integer,                 intent(out) :: i_elm_axis   !< Bezier element, axis is located in
integer,                 intent(out) :: ifail        !< Error code

! --- Local variables
real*8  :: grad_psi, ps_x, ps_y, ps_s, ps_t, xjac
real*8  :: psi_min, psi_max, grad_psi_min
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
integer :: ij_axis(2), i, iv, ms, mt, kf, kv
real*8  :: x(2), s, t, xerr, ferr, rs_tolerance
logical :: early_exit
parameter (rs_tolerance = 1.d-8)

if (my_id .eq. 0) then
  write(*,*) '*********************************'
  write(*,*) '*     find_axis                 *'
  write(*,*) '*********************************'
endif

i_elm_axis = 1
ij_axis    = 1
psi_axis   = 1.d20
grad_psi_min = 1.d20

do i=1,element_list%n_elements

  do ms = 1, 4           ! 4 Gaussian points
    do mt = 1, 4         ! 4 Gaussian points

      ps_s = 0.d0
      ps_t = 0.d0
      R_s  = 0.d0
      Z_s  = 0.d0
      R_t  = 0.d0
      Z_t  = 0.d0
      R    = 0.d0
      Z    = 0.d0

      do kf = 1, 4       ! 4 basis functions
        do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          ps_s = ps_s + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          ps_t = ps_t + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

          R   = R   + node_list%node(iv)%x(kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
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

      if (grad_psi .lt. grad_psi_min) then
        if (     ((tokamak_device(1:4) .ne. 'MAST') .and. (abs(Z) .lt. 0.1 * R_geo)) &
            .or. ((tokamak_device(1:4) .eq. 'MAST') .and. ((abs(Z) .lt. 0.2d0) .and. (R .lt. 1.d0))) ) then
          grad_psi_min = grad_psi
          i_elm_axis = i
          ij_axis(1) = ms;         ij_axis(2)  = mt
        endif
      endif

    enddo
  enddo

enddo

s=Xgauss(ij_axis(1)) ; t=Xgauss(ij_axis(2))

call mnewtax(node_list,element_list,i_elm_axis,s,t,xerr,ferr,ifail)

if ((ifail .ne. 0 ).and.(my_id .eq. 0)) write(*,*) ' MNEWTAX : ifail = ',ifail

call interp(node_list,element_list,i_elm_axis,1,1,s,t,psi_axis,P_s,P_t,P_st,P_ss,P_tt)

call interp_RZ(node_list,element_list,i_elm_axis,s,t,R_axis,R_s,R_t,R_st,R_ss,R_tt,Z_axis,Z_s,Z_t,Z_st,Z_ss,Z_tt)

s_axis = s
t_axis = t

if (my_id .eq. 0) write(*,'(A,i6,4f14.8)') ' magnetic axis : ',i_elm_axis,R_axis,Z_axis,psi_axis

return
end subroutine find_axis

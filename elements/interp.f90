!> subroutine calculates the interpolation within one element (i_elm) for a given position
!! (s,t) in the local coordinates
recursive subroutine interp(node_list, element_list, i_elm, i_var, i_harm, s, t, P, P_s, P_t, P_st,&
  P_ss, P_tt)

use data_structure

implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: i_var
integer,                  intent(in)  :: i_harm
real*8,                   intent(in)  :: s
real*8,                   intent(in)  :: t
real*8,                   intent(out) :: P
real*8,                   intent(out) :: P_s
real*8,                   intent(out) :: P_t
real*8,                   intent(out) :: P_st
real*8,                   intent(out) :: P_ss
real*8,                   intent(out) :: P_tt

! --- Local variables
real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
integer :: kv, iv, kf 

call basisfunctions2(s,t,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),G_ss(1:4,1:4),G_tt(1:4,1:4))

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_st = 0.d0; P_ss = 0.d0; P_tt = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices

  iv = element_list%element(i_elm)%vertex(kv)  ! the node number

  do kf = 1, n_order+1       ! 4 basis functions

    if (i_var == 456) then
#ifdef fullmhd
      P    = P    + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
      P_s  = P_s  + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
      P_t  = P_t  + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
      P_st = P_st + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
      P_ss = P_ss + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
      P_tt = P_tt + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)
#endif
    else
      P    = P    + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
      P_s  = P_s  + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
      P_t  = P_t  + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
      P_st = P_st + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
      P_ss = P_ss + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
      P_tt = P_tt + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)
    endif

  end do

end do

return
end subroutine interp

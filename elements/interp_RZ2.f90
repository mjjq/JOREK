!> Same as interp_RZ, but no second derivatives.
subroutine interp_RZ2(node_list,element_list,i_elm,s,t,R,R_s,R_t,Z,Z_s,Z_t)

use data_structure
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
real*8,                   intent(in)  :: s,t
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

! --- Local variables
real*8  :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
real*8  :: xx1, xx2, ss
integer :: kv, iv, kf

call basisfunctions2(s,t,G,G_s,G_t,G_st,G_ss,G_tt)

R = 0.d0; R_s = 0.d0; R_t = 0.d0
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices

  iv = element_list%element(i_elm)%vertex(kv)  ! the node number

  do kf = 1, n_order+1       ! 4 basis functions
    
    xx1 = node_list%node(iv)%x(kf,1)
    xx2 = node_list%node(iv)%x(kf,2)
    ss  = element_list%element(i_elm)%size(kv,kf)
    
    R    = R    + xx1 * ss * G(kv,kf)
    R_s  = R_s  + xx1 * ss * G_s(kv,kf)
    R_t  = R_t  + xx1 * ss * G_t(kv,kf)

    Z    = Z    + xx2 * ss * G(kv,kf)
    Z_s  = Z_s  + xx2 * ss * G_s(kv,kf)
    Z_t  = Z_t  + xx2 * ss * G_t(kv,kf)

  end do

end do

return
end subroutine interp_RZ2

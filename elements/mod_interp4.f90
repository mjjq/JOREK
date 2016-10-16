module mod_interp4
contains
!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp4(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P)
use data_structure
use phys_module, only: mode
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v)

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), ss
integer :: kv, iv, kf, m, i, i_harm, i_tor

call basisfunctions3(s,t,H,H_s,H_t)

P = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  do kf = 1, n_order+1       ! 4 basis functions
    ss  = element_list%element(i_elm)%size(kv,kf)
    do i = 1, n_v
      P(i)    = P(i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H(kv,kf)
      do i_tor = 1, (n_tor-1)/2
        i_harm = 2*i_tor
        P(i)    = P(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i))   * ss * H(kv,kf)   * cos(mode(i_harm)*phi)
	P(i)    = P(i)   + node_list%node(iv)%values(i_harm+1,kf,i_v(i)) * ss * H(kv,kf)   * sin(mode(i_harm+1)*phi)
      end do
    end do
  end do
end do
return
end subroutine interp4
end module mod_interp4

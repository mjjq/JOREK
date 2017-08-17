module mod_interp_PRZ
contains
!> This subroutine interpolates some variables at positions within one element.
!> Assumes a very simple mode numbering! Namely that every nth mode is a simple multiple of mode 2
pure subroutine interp_PRZ_vec(node_list, element_list, i_elm, i_v, n_v, n_p, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use phys_module, only : mode
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),      intent(in)  :: node_list
type (type_element_list),   intent(in)  :: element_list
integer,                    intent(in)  :: i_elm
integer,                    intent(in)  :: n_v, i_v(n_v), n_p
real*8, dimension(n_p),     intent(in)  :: s, t, phi
real*8, dimension(n_p,n_v), intent(out) :: P, P_s, P_t, P_phi
real*8, dimension(n_p),     intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0

! --- Local variables
real*8  :: H(n_p,4,4), H_s(n_p,4,4), H_t(n_p,4,4), xx1(n_p), xx2(n_p), ss(n_p)
real*8  :: e1r(n_p), e1i(n_p), enr(n_p), eni(n_p), enrtmp(n_p) !< complex and real parts
real*8  :: coss(n_p,n_mode), sins(n_p, n_mode)
integer :: kv, iv, kf, i, i_harm, i_mode, mymodes(n_mode)

call basisfunctions5(n_p,s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0
R = 0.d0; R_s = 0.d0; R_t = 0.d0;
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0;

!XXXX
mymodes = mode(3:n_tor:2)
e1r = cos(mode(2)*phi)
e1i = sin(mode(2)*phi)
enr = 1.d0
eni = 0.d0
do i=1,n_mode
  enrtmp = enr*e1r - eni*e1i
  eni = eni*e1r + enr*e1i
  enr = enrtmp
  coss(:,i) = enr
  sins(:,i) = eni
end do

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  do kf = 1, n_order+1       ! 4 basis functions

    xx1 = node_list%node(iv)%x(kf,1)
    xx2 = node_list%node(iv)%x(kf,2)
    ss  = element_list%element(i_elm)%size(kv,kf)

    R    = R    + xx1 * ss * H(:,kv,kf)
    R_s  = R_s  + xx1 * ss * H_s(:,kv,kf)
    R_t  = R_t  + xx1 * ss * H_t(:,kv,kf)

    Z    = Z    + xx2 * ss * H(:,kv,kf)
    Z_s  = Z_s  + xx2 * ss * H_s(:,kv,kf)
    Z_t  = Z_t  + xx2 * ss * H_t(:,kv,kf)

    do i = 1, n_v
      P(:,i)     = P(:,i)     + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H(:,kv,kf)
      P_s(:,i)   = P_s(:,i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H_s(:,kv,kf)
      P_t(:,i)   = P_t(:,i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H_t(:,kv,kf)
      do i_mode = 1, n_mode

        i_harm = i_mode*2
        P(:,i)     = P(:,i)     + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(:,kv,kf)   * coss(:,i_mode)
        P_s(:,i)   = P_s(:,i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_s(:,kv,kf) * coss(:,i_mode)
        P_t(:,i)   = P_t(:,i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_t(:,kv,kf) * coss(:,i_mode)
        P_phi(:,i) = P_phi(:,i) + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(:,kv,kf)   * sins(:,i_mode)*(-mymodes(i_mode))

        i_harm = i_mode*2+1
        P(:,i)     = P(:,i)     + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(:,kv,kf)   * sins(:,i_mode)
        P_s(:,i)   = P_s(:,i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_s(:,kv,kf) * sins(:,i_mode)
        P_t(:,i)   = P_t(:,i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_t(:,kv,kf) * sins(:,i_mode)
        P_phi(:,i) = P_phi(:,i) + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(:,kv,kf)   * coss(:,i_mode)*mymodes(i_mode)
      enddo
    enddo
  enddo
enddo
end subroutine interp_PRZ_vec

!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use phys_module, only : mode
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
real*8,                   intent(out) :: P_phi(n_v)

integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), xx1, xx2, ss, coss(n_mode), sins(n_mode)
integer :: kv, iv, kf, i, i_harm, i_mode, mymodes(n_mode)
real*8  :: e1r, e1i, enr, eni, enrtmp

call basisfunctions3(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0
R = 0.d0; R_s = 0.d0; R_t = 0.d0;
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0;

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
mymodes = mode(3:n_tor:2)
e1r = cos(mode(2)*phi)
e1i = sin(mode(2)*phi)
enr = 1.d0
eni = 0.d0
do i=1,n_mode
  enrtmp = enr*e1r - eni*e1i
  eni = eni*e1r + enr*e1i
  enr = enrtmp
  coss(i) = enr
  sins(i) = eni
end do
!coss = cos(mymodes*phi)
!sins = sin(mymodes*phi)

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  do kf = 1, n_order+1       ! 4 basis functions

    xx1 = node_list%node(iv)%x(kf,1)
    xx2 = node_list%node(iv)%x(kf,2)
    ss  = element_list%element(i_elm)%size(kv,kf)

    R    = R    + xx1 * ss * H(kv,kf)
    R_s  = R_s  + xx1 * ss * H_s(kv,kf)
    R_t  = R_t  + xx1 * ss * H_t(kv,kf)

    Z    = Z    + xx2 * ss * H(kv,kf)
    Z_s  = Z_s  + xx2 * ss * H_s(kv,kf)
    Z_t  = Z_t  + xx2 * ss * H_t(kv,kf)

    do i = 1, n_v
      P(i)     = P(i)     + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H(kv,kf)
      P_s(i)   = P_s(i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H_s(kv,kf)
      P_t(i)   = P_t(i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H_t(kv,kf)
      do i_mode = 1, n_mode
        i_harm = i_mode*2
        P(i)     = P(i)     + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * coss(i_mode)
        P_s(i)   = P_s(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_s(kv,kf) * coss(i_mode)
        P_t(i)   = P_t(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_t(kv,kf) * coss(i_mode)
        P_phi(i) = P_phi(i) + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * sins(i_mode)*(-mymodes(i_mode))

        i_harm = i_mode*2+1
        P(i)     = P(i)     + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * sins(i_mode)
        P_s(i)   = P_s(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_s(kv,kf) * sins(i_mode)
        P_t(i)   = P_t(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H_t(kv,kf) * sins(i_mode)
        P_phi(i) = P_phi(i) + node_list%node(iv)%values(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * coss(i_mode)*mymodes(i_mode)
      enddo
    enddo
  enddo
enddo
end subroutine interp_PRZ


!> This subroutine interpolates the change (u_n-u_(n-1)) in some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ_delta(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use phys_module, only : mode
use mod_basisfunctions
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
real*8,                   intent(out) :: P_phi(n_v)

integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), xx1, xx2, ss, coss(n_mode), sins(n_mode)
integer :: kv, iv, kf, i, i_harm, i_mode, mymodes(n_mode)

call basisfunctions3(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0
R = 0.d0; R_s = 0.d0; R_t = 0.d0;
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0;

mymodes = mode(3:n_tor:2)
coss = cos(mymodes*phi)
sins = sin(mymodes*phi)

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  do kf = 1, n_order+1       ! 4 basis functions

    xx1 = node_list%node(iv)%x(kf,1)
    xx2 = node_list%node(iv)%x(kf,2)
    ss  = element_list%element(i_elm)%size(kv,kf)

    R    = R    + xx1 * ss * H(kv,kf)
    R_s  = R_s  + xx1 * ss * H_s(kv,kf)
    R_t  = R_t  + xx1 * ss * H_t(kv,kf)

    Z    = Z    + xx2 * ss * H(kv,kf)
    Z_s  = Z_s  + xx2 * ss * H_s(kv,kf)
    Z_t  = Z_t  + xx2 * ss * H_t(kv,kf)

    do i = 1, n_v
      P(i)     = P(i)     + node_list%node(iv)%deltas(1,kf,i_v(i)) * ss * H(kv,kf)
      P_s(i)   = P_s(i)   + node_list%node(iv)%deltas(1,kf,i_v(i)) * ss * H_s(kv,kf)
      P_t(i)   = P_t(i)   + node_list%node(iv)%deltas(1,kf,i_v(i)) * ss * H_t(kv,kf)
      do i_mode = 1, n_mode
        i_harm = i_mode*2
        P(i)     = P(i)     + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * coss(i_mode)
        P_s(i)   = P_s(i)   + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H_s(kv,kf) * coss(i_mode)
        P_t(i)   = P_t(i)   + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H_t(kv,kf) * coss(i_mode)
        P_phi(i) = P_phi(i) + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * sins(i_mode)*(-mymodes(i_mode))

        i_harm = i_mode*2+1
        P(i)     = P(i)     + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * sins(i_mode)
        P_s(i)   = P_s(i)   + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H_s(kv,kf) * sins(i_mode)
        P_t(i)   = P_t(i)   + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H_t(kv,kf) * sins(i_mode)
        P_phi(i) = P_phi(i) + node_list%node(iv)%deltas(i_harm,kf,i_v(i)) * ss * H(kv,kf)   * coss(i_mode)*mymodes(i_mode)
      enddo
    enddo
  enddo
enddo
end subroutine interp_PRZ_delta
end module mod_interp_PRZ

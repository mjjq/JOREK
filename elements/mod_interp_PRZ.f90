module mod_interp_PRZ
contains
!> This subroutine interpolates some variables at positions within one element.
!> Assumes a very simple mode numbering! Namely that every nth mode is a simple multiple of mode 2
!>
!> This is slightly faster for many particles since it only needs to precalculate
!> the values once. Save 10-20% for > 10 particles in one call.
pure subroutine interp_PRZ_vec(node_list, element_list, i_elm, i_v, n_v, n_p, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use mod_parameters, only: n_period
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
real*8  :: H(4,4), H_s(4,4), H_t(4,4), HZ(n_tor,n_p), dHZ(n_tor,n_p)
real*8  :: e1r(n_p), e1i(n_p), enr(n_p), eni(n_p), enrtmp(n_p) !< complex and real parts
integer :: kv, iv, kf, i, j, i_harm, i_mode, i_tor, i_p
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
e1r = cos(n_period*phi)
e1i = sin(n_period*phi)
enr = 1.d0
eni = 0.d0
HZ(1,:) = 1.d0
dHZ(1,:) = 0.d0
do i=1,n_mode
  enrtmp = enr
  enr = enr*e1r - eni*e1i
  eni = eni*e1r + enrtmp*e1i
  HZ(2*i,:)   = enr
  HZ(2*i+1,:) = eni
end do
do i=1,n_mode
  dHZ(2*i,:)   = HZ(2*i+1,:)*(-n_period*i)
  dHZ(2*i+1,:) = HZ(2*i,:)*(n_period*i)
end do

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  do i = 1, n_v
    do kf=1,n_order+1
      values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
    end do
  end do
  xR(:,kv) = node_list%node(iv)%x(:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(:,2) * sizes(:)
end do

do i_p = 1,n_p
  call basisfunctions6(s(i_p),t(i_p),H,H_s,H_t)
  R(i_p)   = sum(xR*H)
  R_s(i_p) = sum(xR*H_s)
  R_t(i_p) = sum(xR*H_t)
  Z(i_p)   = sum(xZ*H)
  Z_s(i_p) = sum(xZ*H_s)
  Z_t(i_p) = sum(xZ*H_t)

  do kv = 1,n_vertex_max
    do i = 1, n_v
      do kf = 1, n_order+1
        v = sum(values(1:n_tor,kf,i,kv)*HZ(1:n_tor, i_p))
        P(i_p,i)     = P(i_p,i)     + v * H(kf, kv)
        P_s(i_p,i)   = P_s(i_p,i)   + v * H_s(kf, kv)
        P_t(i_p,i)   = P_t(i_p,i)   + v * H_t(kf, kv)
        vp = sum(values(1:n_tor,kf,i,kv)*dHZ(1:n_tor, i_p))
        P_phi(i_p,i) = P_phi(i_p,i) + vp * H(kf, kv)
      enddo
    enddo
  enddo
enddo
end subroutine interp_PRZ_vec

!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use mod_parameters, only: n_period
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
real*8  :: H(4,4), H_s(4,4), H_t(4,4), HZ(n_tor), dHZ(n_tor)
integer :: kv, iv, kf, i, i_tor
real*8  :: e1r, e1i, enr, eni, enrtmp
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp

!maybe try !DIR$ FORCEINLINE 
call basisfunctions6(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
e1r = cos(n_period*phi)
e1i = sin(n_period*phi)
enr = 1.d0
eni = 0.d0
HZ(1) = 1.d0
dHZ(1) = 0.d0
do i=1,n_mode
  enrtmp = enr
  enr = enr*e1r - eni*e1i
  eni = eni*e1r + enrtmp*e1i
  HZ(2*i)   = enr
  HZ(2*i+1) = eni
  dHZ(2*i) = eni*(-n_period*i)
  dHZ(2*i+1) = enr*(n_period*i)
end do

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  do i = 1, n_v
    do kf=1,n_order+1
      values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
    end do
  end do
  xR(:,kv) = node_list%node(iv)%x(:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(:,2) * sizes(:)
end do

R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)

do kv = 1,n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = sum(values(1:n_tor,kf,i,kv)*HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
      P_s(i)   = P_s(i)   + v * H_s(kf, kv)
      P_t(i)   = P_t(i)   + v * H_t(kf, kv)
      vp = sum(values(1:n_tor,kf,i,kv)*dHZ(1:n_tor))
      P_phi(i) = P_phi(i) + vp * H(kf, kv)
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

subroutine sputtering(i_spec, E_in,angle_in,Y,angle_out,E_out)

!use constants, only : EL_CHG, PI

implicit none

integer :: i_spec
real*8  :: E_in, angle_in
real*8  :: Y, E_out, angle_out

real*8 :: E_0, eps_L, eps_L2, a_L, a_B, s_n_KrC, w_eps, E_frac, Y_0, alfa, alfa_0, cos_alfa

! sputtering formula from: W. Eckstein / Vacuum 82 (2008) 930–934
! Behrisch, Eckstein, Sputtering by Particle Bombardment, Springer, Topics in Applied Physics 110, p142

real*8, dimension(3), parameter :: Lambda = (/ 1.0087d0, 0.3583d0, 0.2879d0 /)      ! H,D,T (1) on W (2)
real*8, dimension(3), parameter :: q      = (/ 0.0075d0, 0.0183d0, 0.0419d0 /)
real*8, dimension(3), parameter :: zmu    = (/ 1.2046d0, 1.4410d0, 1.58202d0 /)
real*8, dimension(3), parameter :: E_th   = (/ 457.42d0, 228.84d0, 153.8842d0/)     ! [eV]
real*8, dimension(3), parameter :: eps    = (/ 9.86986d3,9.92326d3, 9.97718d3 /)

real*8, dimension(3), parameter :: m_1 = (/ 1.d0, 2.d0, 3.d0 /)      ! projectile atomic mass
real*8, dimension(3), parameter :: Z_1 = 1.d0                        !            atomic number
real*8,               parameter :: m_2 = 183.84d0                    ! target     atomic mass (Tungsten)
real*8,               parameter :: Z_2 = 74.d0                       !            atomic number

real*8, parameter :: E_sb       = 8.68
real*8, parameter :: E_sb_gamma = 202.85
real*8, parameter :: e2_eVm     = 14.4d-10   ! [eV m]

real*8, dimension(3), parameter :: f = (/ 1.3708d0, 1.1544d0, 1.1499d0 /)
real*8, dimension(3), parameter :: b = (/ 0.4824d0, 0.1901d0, 0.1573d0 /)
real*8, dimension(3), parameter :: c = (/ 1.0067d0, 1.0824d0, 1.1050d0 /)
real*8, parameter :: E_sp = 1.0 !eV

real*8, parameter :: PI = 3.1415926535897

E_0   = E_in
eps_L = E_0 / eps(i_spec)

!a_B   = 0.0529177d-9                                           ! Bohr radius [m]
!a_L   = 0.885341377 * a_B / sqrt( Z_1(i_spec)**0.66666 + Z_2**0.66666) ! Lindhard screening length
!eps_L2 = E_0 * m_2 / (m_1(i_spec)+m_2) * a_L /(Z_1(i_spec) * Z_2 * e2_eVm)

w_eps = eps_L + 0.1728d0 * sqrt(eps_L) + 0.008d0 * eps_L**0.1504d0

s_n_KrC = 0.5d0 * Log(1.d0 + 1.2288*eps_L) / w_eps

E_frac = max(E_0/E_th(i_spec) - 1.d0,0.d0)**zmu(i_spec)

Y_0 = q(i_spec) * s_n_KrC * E_frac / ( lambda(i_spec) / w_eps + E_frac )

alfa   = angle_in
alfa_0 = max(PI/2.d0, PI - acos(sqrt(1.d0/(1.d0 + E_0/E_sp))))

cos_alfa = cos((alfa/alfa_0 * PI/2.d0)**c(i_spec))

Y = Y_0 / cos_alfa**f(i_spec) * exp( b(i_spec) * (1.d0 - 1.d0 / cos_alfa))

E_out = E_sb /2.d0  ! cosine**2(theta) distribution

write(*,*) Y,Y_0


return
end

!program try
!real*8 :: E_in

!E_in = 100.
!angle_in = 3.1415926 / 3.

!call  sputtering(2,E_in,angle_in,Y,angle_out,E_out)

!end

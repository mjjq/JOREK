!> Contains the Bezier and Fourier basis functions on Gaussian points.
!!
!! - Array dimensions of Bezier basis functions and derivatives
!!   (H, H_s, H_t, H_st, H_ss, H_tt):
!!   - Vertex
!!   - Basis function (i.e. p0,u,v,w)
!!   - Gaussian point in s direction
!!   - Gaussian point in t direction
!! - Array dimensions of Fourier basis functions and derivatives
!!   (HZ, HZ_p):
!!   - Toroidal mode index
!!   - Toroidal plane index
!! - Array dimensions of one-dimensional basis functions
!!   (H1, H1_s, H1_ss)
!!   - Vertex
!!   - Basis function
!!   - Gaussian point
!!
!! @see Gauss
module basis_at_gaussian
  
  use mod_parameters
  use gauss
  
  implicit none
  
  save
  
  real*8 :: H(n_vertex_max,n_order+1,n_gauss,n_gauss)    !< Basis functions in poloidal plane
  real*8 :: H_s(n_vertex_max,n_order+1,n_gauss,n_gauss)  !< Basis functions in poloidal plane; derivative with respect to first coordinate
  real*8 :: H_t(n_vertex_max,n_order+1,n_gauss,n_gauss)  !< Basis functions in poloidal plane; derivative with respect to second coordinate
  real*8 :: H_st(n_vertex_max,n_order+1,n_gauss,n_gauss) !< Basis functions in poloidal plane; cross-derivative 
  real*8 :: H_ss(n_vertex_max,n_order+1,n_gauss,n_gauss) !< Basis functions in poloidal plane; second derivative with respect to first coordinate
  real*8 :: H_tt(n_vertex_max,n_order+1,n_gauss,n_gauss) !< Basis functions in poloidal plane; second derivative with respect to second coordinate

  real*8 :: HZ(n_tor,n_plane)          !< Basis functions in toroidal direction
  real*8 :: HZ_p(n_tor,n_plane)        !< Derivative of basis functions in toroidal direction
  real*8 :: HZ_pp(n_tor,n_plane)       !< Second derivative of basis functions in toroidal direction

  real*8 :: H1(2,n_order+1,n_gauss)    !< One dimensional basis functions
  real*8 :: H1_s(2,n_order+1,n_gauss)  !< First derivative of one dimensional basis functions
  real*8 :: H1_ss(2,n_order+1,n_gauss) !< Second derivative of one dimensional basis functions
  
end module basis_at_gaussian

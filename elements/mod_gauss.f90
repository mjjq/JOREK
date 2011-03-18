!> Contains positions (xgauss) and weights (wgauss) of
!! Gaussian points for Gaussian integration.
!!
!! The values are valid for normalised coordinates in the range
!! \f$0 \le S \le 1\f$.
module Gauss
  
 integer, parameter :: n_gauss   = 4                  !< Number of Gaussian points
 integer, parameter :: n_gauss_2 = n_gauss * n_gauss  !< Square of n_gauss
 
 real*8,  parameter :: Xgauss(n_gauss) = (/ 0.0694318442029735d0, 0.3300094782075720d0,            &
   0.6699905217924280d0, 0.9305681557970265d0 /)      !< Positions of Gaussian points

 real*8,  parameter :: Wgauss(n_gauss) = (/ 0.173927422568727d0,  0.326072577431273d0,             &
   0.326072577431273d0,  0.173927422568727d0  /)      !< Weights of Gaussian points

end module Gauss

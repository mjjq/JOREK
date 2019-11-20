!> Contains positions (xgauss) and weights (wgauss) of
!! Gaussian points for Gaussian integration.
!!
!! The values are valid for normalised coordinates in the range
!! \f$0 \le S \le 1\f$.
!!
!! Taken from https://pomax.github.io/bezierinfo/legendre-gauss.html and
!! converted to [0,1], with weights normalized to sum 1
module gauss

#if GAUSS_ORDER == 8
 integer, parameter :: n_gauss   = 8                  !< Number of Gaussian points
 
 real*8,  parameter :: Xgauss(n_gauss) = &
     [0.019855071751232d0, 0.101666761293187d0, 0.237233795041836d0, 0.408282678752176d0, 0.591717321247825d0, &
     0.762766204958165d0, 0.898333238706813d0, 0.980144928248768d0]

 real*8,  parameter :: Wgauss(n_gauss) = &
     [0.050614268145188d0, 0.111190517226687d0, 0.156853322938944d0, 0.181341891689181d0, 0.181341891689181d0, &
     0.156853322938944d0, 0.111190517226687d0, 0.050614268145188d0]
#else
  integer, parameter :: n_gauss   = 4                  !< Number of Gaussian points
  real*8,  parameter :: Xgauss(n_gauss) = &
    [0.0694318442029735d0, 0.3300094782075720d0, 0.6699905217924280d0, 0.9305681557970265d0]      !< Positions of Gaussian points

  real*8,  parameter :: Wgauss(n_gauss) = &
    [0.173927422568727d0,  0.326072577431273d0, 0.326072577431273d0,  0.173927422568727d0]      !< Weights of Gaussian points
#endif

  integer, parameter :: n_gauss_2 = n_gauss * n_gauss  !< Square of n_gauss
end module gauss

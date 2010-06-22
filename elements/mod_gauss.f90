module Gauss
!-----------------------------------------------------------
! module contains the Gaussian points (Xgauss) and weights
! (Wgauss) in the normalised coordinates ( 0. <= S <= 1.)
!-----------------------------------------------------------
 integer, parameter :: n_gauss = 4, n_gauss_2 = n_gauss * n_gauss
 real*8 :: Xgauss(n_gauss), Wgauss(n_gauss)

 data XGauss /0.0694318442029735,0.3300094782075720,0.6699905217924280,0.9305681557970265 /
 data WGauss /0.173927422568727, 0.326072577431273, 0.326072577431273 ,0.173927422568727/

endmodule

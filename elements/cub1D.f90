subroutine cub1D(X1,X1S,X2,X2S,S,X,XS)
!-----------------------------------------------------------------------
! CUBIC HERMITE INTERPOLATION IN ONE DIMENSION
!-----------------------------------------------------------------------
implicit none
real*8   :: X1,X1S,X2,X2S,S,X,XS
real*8   :: H0M,H0P,H1M,H1P,H0MS,H0PS,H1MS,H1PS

H0M  =  (S-1.)**2 *(S+2.) * 0.25
H0MS =  (S-1.)*(S+2.)/2. + (S-1.)**2 * 0.25
H0P  = -(S+1.)**2 *(S-2.) * 0.25
H0PS = -(S+1.)*(S-2.)/2. - (S+1.)**2 * 0.25
H1M  =  (S-1.)**2 *(S+1.) * 0.25
H1MS =  (S-1.)*(S+1.)/2. + (S-1.)**2 * 0.25
H1P  =  (S+1.)**2 *(S-1.) * 0.25
H1PS =  (S+1.)*(S-1.)/2. + (S+1.)**2 * 0.25

X  = X1*H0M  + X1S*H1M +  X2*H0P  + X2S*H1P
XS = X1*H0MS + X1S*H1MS + X2*H0PS + X2S*H1PS

return
end

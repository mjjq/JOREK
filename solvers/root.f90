function root(A,B,C,D,SGN)
!---------------------------------------------------------------------
! THIS FUNCTION GIVES BETTER ROOTS OF QUADRATICS BY AVOIDING
! CANCELLATION OF SMALLER ROOT
!---------------------------------------------------------------------
implicit none
real*8 :: root,a, b, c, d, sgn

if(((B .EQ. 0.D0) .and. (D .EQ. 0.D0)) .or. (A .eq. 0.D0)) then
 root = 1.d20 ! ill defined
 return
endif
 
if (B*SGN .GE. 0.d0) then
  root = -2.d0*C/(B+SGN*SQRT(D))
else
  ROOT = (-B + SGN*SQRT(D)) / (2.d0 * A)
endif
return
end

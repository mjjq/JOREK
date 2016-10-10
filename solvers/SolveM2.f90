subroutine SolveM2(a,b,c,d,e,f,x,y)
!-----------------------------------------------------------------------
! solves a 2x2 set of equations
!-----------------------------------------------------------------------
implicit none
real*8 :: a,b,c,d,e,f,x,y, det

det = a*d - b*c

if (det .ne. 0.) then
  x = (   d * e - b * f ) / det
  y = ( - c * e + a * f ) / det
else
  x = 0.
  y = 0.
endif

return
end
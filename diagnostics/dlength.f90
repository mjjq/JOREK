real*8 function dlength(x1,x2)
use mod_parameters
real*8 :: total, x1(*), x2(*)
integer :: i
total = 0.
do i=1,n_dim
 total = total + (x1(i)-x2(i))**2
enddo
dlength = sqrt(total)

return
end

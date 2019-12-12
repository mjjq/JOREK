subroutine plot_coils(frame)

use vacuum_response
use phys_module, only: write_ps

implicit none

! --- Routine parameters
logical, intent(in) :: frame

! --- Local variables
real*8  :: rp(5),zp(5), r_min, r_max, z_min, z_max
integer :: i

if ( .not. write_ps  ) then
  write(*,*) 'Jorek2postscript deactivated. Skipping plot_coils'
  return
endif
write(*,*) '*************************************'
write(*,*) '*         plot coils                *'
write(*,*) '*************************************'

R_max = maxval(R_coils+dR_coils)
R_min = minval(R_coils-dR_coils)
Z_max = maxval(Z_coils+dZ_coils)
Z_min = minval(Z_coils-dZ_coils)

R_min = max(0.d0,R_min)

if (frame) CALL NFRAME(21,11,1,R_min,R_max,Z_min,Z_max,' ',1,'R [m]',5,'Z [m]',5)

do i=1,n_coils

  rp(1) = R_coils(i) - dR_coils(i)/2.
  rp(2) = R_coils(i) + dR_coils(i)/2.
  rp(3) = R_coils(i) + dR_coils(i)/2.
  rp(4) = R_coils(i) - dR_coils(i)/2.
  rp(5) = rp(1)

  zp(1) = Z_coils(i) - dZ_coils(i)/2.
  zp(2) = Z_coils(i) - dZ_coils(i)/2.
  zp(3) = Z_coils(i) + dZ_coils(i)/2.
  zp(4) = Z_coils(i) + dZ_coils(i)/2.
  zp(5) = zp(1)
  
  write(51,*) ' 2. setlinewidth '
  call lplot6(21,11,rp,zp,-5,' ')

enddo
write(51,*) ' .5 setlinewidth '

return
end

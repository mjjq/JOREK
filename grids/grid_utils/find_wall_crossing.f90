subroutine find_wall_crossing(R_wall,Z_wall,n_wall,Rp,Zp,tht_p,Rw,Zw,Tw)
  !-----------------------------------------------------------------------
  ! subroutine to find the crossing of the line given by (Rp,Zp,tht_p)
  ! with one of the wall segments (given by straight lines, R_wall,Z_wall)
  ! The result is the position of the crossing (Rw,Zw)
  ! (only the first found crossing is returned)
  !-----------------------------------------------------------------------
  
  implicit none
  
  ! --- input/output variables
  real*8,  intent(in)  :: R_wall(*), Z_wall(*)
  real*8               :: Rp,        Zp,       tht_p
  real*8,  intent(out) :: Rw,        Zw,       Tw
  integer, intent(in)  :: n_wall
  
  ! --- local variables
  real*8  :: tan_p, tan12, PI, TWOPI, HALFPI, area
  integer :: k, my_id
  logical :: found
  
  tan_p  = tan(tht_p)
  HALFPI = asin(1.d0)
  PI     = 2.d0 * HALFPI
  TWOPI  = 2.d0 * PI
  my_id  = 0
  
  Rw = 0.d0; ZW = 0.d0
  
  tht_p = mod(tht_p,TWOPI)
  if (tht_p .lt. -HALFPI) tht_p = tht_p + TWOPI
  if (tht_p .gt. 1.5*PI)  tht_p = tht_p - TWOPI
  
  found = .false.

!---- determine if wall points are clockwise or anti-clockwise
  area = 0.d0           ! assumes a closed curve (first point equals last point)
  do k=1,n_wall-1
    area = area + (R_wall(k+1)-R_wall(k))*(Z_wall(k+1)+Z_wall(k))
  enddo
  
  do k=1,n_wall-1
  
    if (R_wall(k) .ne. R_wall(k+1)) then

      tan12 = (Z_wall(k+1) - Z_wall(k))/(R_wall(k+1) - R_wall(k))

      if (Z_wall(k) .ne. Z_wall(k+1)) then

        Rw    = (Z_wall(k) - R_wall(k)*tan12 - Zp + Rp*tan_p ) / (tan_p - tan12)

        Zw    = Zp + tan_p  * (Rw - Rp)

      elseif (tan_p .ne. 0.d0) then

        Zw = Z_wall(k)

        Rw = Rp + (Zw - Zp)/ tan_p

      else

        Rw = 1.d20  ! no solution
        Zw = 1.d20

      endif

    else

      Rw = R_wall(k)

      Zw = Zp + tan_p * (Rw - Rp)

    endif
  
    if (((Zw-Z_wall(k))*(Zw-Z_wall(k+1)) .le. 0.d0) .and. (((Rw-R_wall(k))*(Rw-R_wall(k+1)) .le. 0.d0))) then

      if (area .lt. 0) then
        Tw = atan2(Z_wall(k+1)- Z_wall(k),  R_wall(k+1)- R_wall(k))
      else
        Tw = atan2(Z_wall(k)  - Z_wall(k+1),R_wall(k)  - R_wall(k+1))
      endif
   
      if (((tht_p .gt. -HALFPI) .and. (tht_p .le. HALFPI))    .and. (Rw .ge. Rp)) then

        found = .true.

      elseif (((tht_p .gt. HALFPI) .and. (tht_p .lt. 1.5*PI)) .and. (Rw .le. Rp)) then

        found = .true.

      endif
    endif
  
    if (found) exit
  
  enddo
  
!  if (.not. found) write(*,'(A,6f16.8)') ' CROSSING NOT FOUND : ',Rp,Zp,tht_p,Rw,Zw,Tw
  
  return
end subroutine find_wall_crossing


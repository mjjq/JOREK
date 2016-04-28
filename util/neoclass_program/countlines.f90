module countlines_mod
implicit none
contains
function countlines(fd)
  integer countlines
  integer, intent(in) :: fd
  
  integer i

  i=1
  DO
     READ(fd,'()', END=100)
     i=i+1
  ENDDO
100 REWIND(fd)
  countlines = i-1
end function countlines
end module countlines_mod

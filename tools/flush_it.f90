subroutine flush_it(unit)
  
  implicit none
  
  integer, intent(in) :: unit
  
#ifdef IBM_MACHINE
  call flush_(unit)
#else
  call flush(unit)
#endif
  
end subroutine flush_it

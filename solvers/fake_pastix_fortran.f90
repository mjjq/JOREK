!> A fake pastix_fortran subroutine to prevent linking problems
subroutine fake_pastix_fortran
  write(*,*) "ERROR"
  call exit(100)
end subroutine fake_pastix_fortran

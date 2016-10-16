!> Print some information from version.h
subroutine print_version
#include "version.h"
  111 format(2x,a,': ',a)
  write(*,*) ' ', trim(adjustl(RCS_VERSION))
  write(*,111) 'compile_time        ', trim(adjustl(compile_time))
  write(*,111) 'compile_user        ', trim(adjustl(compile_user))
  write(*,111) 'compile_machine     ', trim(adjustl(compile_machine))
  write(*,111) 'compile_dir         ', trim(adjustl(compile_dir))
  write(*,111) 'compile_command     ', trim(adjustl(compile_command))
  write(*,111) 'compile_flags       ', trim(adjustl(compile_flags))
  write(*,111) 'compile_includes    ', trim(adjustl(compile_includes))
  write(*,111) 'compile_defines     ', trim(adjustl(compile_defines))
  write(*,111) 'compile_libs        ', trim(adjustl(compile_libs))
  write(*,111) 'compile_modules     ', trim(adjustl(compile_modules))
end subroutine print_version

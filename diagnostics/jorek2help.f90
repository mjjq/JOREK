!> Command line function, print help or hardcoded parameter values to screen
!!
!! NOTE: To remove this function a -DNO_HELP flag should be added to the
!! compiler flags
!!
!! When one or more commandline argument is given this function will loop
!! through those options and will respond to each one till one is recognized.
!! The following argument are supported:
!!
!!   -h or --help  : print some basic info about JOREK_2
!!   -p or --param : print hardcoded parameters
!!
!! unknown argument give a short message that the argument is invalit. When
!! none of the given arguments is recognized a sugestion is given for valid
!! arguments.
!!
!! writer: Mark Verbeek
!! data  : 11-02-2016
subroutine jorek2help()
  use parameters
  use mod_log_params

!#include "r3_info.h"
#include "version.h"
  
  implicit none

#ifndef NO_HELP

  integer           :: narg,cptArg  !> for commandline arguments
  character(len=20) :: ArgName      !> Argument name

  narg = command_argument_count() !> get number of commandline agruments

  !> when argument is give print info
  111 format('*',2x,a,': ',a)
  200 format(79('*'))
    
  if(narg>0)then
    do cptArg=1,narg
      call get_command_argument(cptArg,ArgName)
      select case(adjustl(ArgName))
        case("--about","-a")
          write(*,200)
          write(*,'(A)') '*                    3D Reduced MHD : JOREK_2.0                               *'
          write(*,200)
          write(*,'(A)') '* Solves the (reduced) MHD equations in 3D toroidal geometry                  *'
          write(*,'(A)') '*                                                                             *'
          write(*,'(A)') '* - solvers implemented:                                                      *'
          write(*,'(A)') '*   - MUMPS                                                                   *'
          write(*,'(A)') '*   - PastiX                                                                  *'
          write(*,'(A)') '*   - GMRES (+MUMPS or PastiX preconditioner)                                 *'
          write(*,'(A)') '*                                                                             *'
          write(*,'(A)') '* - required libraries :                                                      *'
          write(*,'(A)') '*   - MPI                                                                     *'
          write(*,'(A)') '*   - MUMPS                                                                   *'
          write(*,'(A)') '*   - PastiX                                                                  *'
          write(*,'(A)') '*   - SCOTCH (metis)                                                          *'
          write(*,'(A)') '*   - FFTW                                                                    *'
          write(*,'(A)') '*   - SCALAPACK (BLACS)                                                       *'
          write(*,'(A)') '*   - LAPACK, BLAS                                                            *'
          write(*,'(A)') '*   - PPPLIB                                                                  *'
          write(*,'(A)') '*                                                                             *'
          write(*,'(A)') '* Original author : Guido Huysmans (Euratom / CEA Association)                *'
          write(*,'(A)') '*          authors: ...                                                       *'
          write(*,'(A)') '* start date: 18-7-2008                                                       *'
          write(*,'(A)') '*                                                                             *'
          write(*,200)
          stop
        case("--help","-h")
          write(*,200)
          write(*,'(A)') '* List of command line argument options                                       *'
          write(*,'(A)') '*   -a or --about for info about JOREK in general                             *'
          write(*,'(A)') '*   -p or --param for JOREK hardcoded parameters                              *'
          write(*,'(A)') '*   -v or --version for JOREK compile veriosn                                 *'
          write(*,200)
          stop
        case("--version","-v")
          write(*,'(A)') '* --- Compile Info                                                            *'
          write(*,*)   '* RCS revision        : ', RCS_VERSION
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
          stop
        case("--param","-p")
          call log_parameters(0, .true.)
          stop
        case default
          write(*,'(A, A, A)') 'Option ', adjustl(ArgName), 'unknown'
      end select
    end do
    write(*,'(A)') 'try -h or -help for program info'
    stop
  end if

#endif

return
end subroutine jorek2help

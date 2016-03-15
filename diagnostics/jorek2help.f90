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

  implicit none

#ifndef NO_HELP

  integer           :: narg,cptArg  !> for commandline arguments
  character(len=20) :: ArgName      !> Argument name

  narg = command_argument_count() !> get number of commandline agruments

  !> when argument is give print info
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

        case("--param","-p")
          call log_parameters(0, .true.)
          stop
        case default
          write(*,'(A, A, A)') 'Option ', adjustl(ArgName), 'unknown'
      end select
    end do
    write(*,'(A)') 'try -h or -help for program info or -p for parameters list'
    stop
  end if

#endif

return
end subroutine jorek2help

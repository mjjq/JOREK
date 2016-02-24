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
  
  integer           ::narg,cptArg   !> for commandline arguments
  integer           :: i, j, n_rows !> do loop index
  character(len=20) :: ArgName      !> Argument name
  
  narg = command_argument_count() !> get number of commandline agruments
  
  !> when argument is give print info
  200 format(79('*'))
  if(narg>0)then    
    do cptArg=1,narg
      call get_command_argument(cptArg,ArgName)
      select case(adjustl(ArgName))
        case("--help","-h")
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
      case("--param","-p")
        112 format(A, i12, 41X, A)
        write(*,200)
        write(*,'(A)') '*                    3D Reduced MHD : JOREK_2.0                               *'
        write(*,200)
        write(*,'(A)') '* Using following parameters:                                                 *'
        write(*,  112) '*      jorek_model    =  ', jorek_model       , '*'
        write(*,  112) '*      n_var          =  ', n_var             , '*'
        write(*,  112) '*      n_dim          =  ', n_dim             , '*'
        write(*,  112) '*      n_order        =  ', n_order           , '*'
        write(*,  112) '*      n_tor          =  ', n_tor             , '*'
        write(*,  112) '*      n_period       =  ', n_period          , '*'
        write(*,  112) '*      n_plane        =  ', n_plane           , '*'
        write(*,  112) '*      n_vertex_max   =  ', n_vertex_max      , '*'
        write(*,  112) '*      n_elements_max =  ', n_elements_max    , '*'
        write(*,  112) '*      n_boundary_max =  ', n_boundary_max    , '*'
        write(*,  112) '*      n_pieces_max   =  ', n_pieces_max      , '*'
        write(*,  112) '*      n_degrees      =  ', n_degrees         , '*'
        write(*,  112) '*      nref_max       =  ', nref_max          , '*'
        write(*,  112) '*      n_ref_list     =  ', n_ref_list        , '*'
        write(*,'(A)') '*                                                                             *' 
        write(*,'(A)') '*      Implemented variables:                                                 *' 
        
        ! determine number of rows needed to show all variable_names
        n_rows = ceiling(n_var/4.0)  
        
        ! The first loop loops through the row needed. The the left eastectics is 
        ! written followed by a loop that print out the variable_name of white space 
        ! depending on it this variable_name exist. The last write is the eastectics 
        ! on the right. 
        do i = 0,n_rows-1
          write(*,'(A)',advance='no') '*      '
          do j = (i*4) + 1, (i*4) + 4
            if ( j .gt. n_var) then
              write(*,'(11x)',advance='no') 
            else
              write(*,'(A11)',advance='no') variable_names(j)
            end if
            if ( j .lt. (i*4 + 4)) then
              write(*,'(7x)',advance='no') 
            end if
          end do 
          write(*,'(A)') '      *'
        end do

        write(*,'(A)') '*                                                                             *' 
        write(*,200) 
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

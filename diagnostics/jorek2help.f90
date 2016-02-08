! Command line function

subroutine jorek2help()
  use parameters
  
  implicit none

  integer::narg,cptArg  !> for commandline arguments
  character(len=20)::ArgName  !> Argument name

  narg = command_argument_count()
  
  !> when argument is give print info of first one and exit normally
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

        113 format(A, A11, A, A11, A, A11, A, A11, A)
        write(*,113) '*      ', variable_names(1), '       ', variable_names(2), '       ', &
                                variable_names(3), '       ', variable_names(4), '      *' 
#if (JOREK_MODEL == 500 || JOREK_MODEL == 555)
        write(*,113) '*      ', variable_names(5), '       ', variable_names(6), '       ', &
                                variable_names(7), '       ', variable_names(8), '      *' 
#else
        write(*,113) '*      ', variable_names(5), '       ', variable_names(6), '       ', &
                                variable_names(7), '       ', ' ',               '      *' 
#endif
        write(*,'(A)') '*                                                                             *' 
        write(*,200) 
        stop
      case default
        write(*,'(A)') 'Option ', adjustl(ArgName), 'unknown, try -h or -help for program info '
        write(*,'(A)') 'or -p for parameters list'
      stop
        end select
    end do
  end if

return  
end subroutine jorek2help

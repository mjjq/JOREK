!**********************************************************************
!* program to extract data from a JOREK2 restart file                 *
!**********************************************************************

program jorek2_diagno
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
use pellet_module
use mpi_mod
use mod_import_restart
use domains
use mod_log_params
use diagnostics, only: axis_is_psi_minimum

use mod_element_rtree

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
integer :: i, in, i_tor, i_spi, resultlength
real*8  :: growth_kin, growth_mag,density,density_in,density_out,pressure,pressure_in,pressure_out
real*8  :: Rplot(2), Zplot(2)
real*8  :: psi_bnd, psi_axis,R_axis,Z_axis,s_axis,t_axis
real*8  :: spi_abl_rate_tot, spi_abl_tot
integer :: ifail, my_id, ierr, i_elm_axis
integer :: required, provided, StatInfo
integer :: rank, comm_size, n_cpu
integer :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS

  interface
    subroutine distribute_vector(my_id,rhs,rhs_dis,again)
      real*8               :: rhs(:), rhs_dis(:)
      integer              :: my_id
      logical              :: again
    end subroutine distribute_vector

    subroutine distribute_harmonics(my_id,my_id_n,n_cpu)
      integer              :: my_id, my_id_n,n_cpu
    end subroutine distribute_harmonics

    subroutine gmres_driver(my_id,my_id_n,i_tor,n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
      integer :: i_tor(:), my_id, my_id_n, MPI_COMM_N, MPI_COMM_MASTER
      integer :: iter_gmres, n_tor
    end subroutine gmres_driver


    subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
      use data_structure
      integer(kind=4),             intent(in)    :: my_id
      integer(kind=4),             intent(in)    :: xcase2
      type (type_node_list),       intent(inout) :: node_list
      type (type_element_list),    intent(inout) :: element_list
      type (type_bnd_node_list)   ,intent(inout) :: bnd_node_list
      type (type_bnd_element_list),intent(inout) :: bnd_elm_list
      logical(kind=4),             intent(in)    :: xpoint2
      logical(kind=4),             intent(in)    :: nice_q
    end subroutine equilibrium

    subroutine set_trap_sigterm() bind(C)
    end subroutine set_trap_sigterm
    logical function sigterm_called() bind(C)
    end function sigterm_called
  end interface


real*8  :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), mindelta, maxdelta
real*8  :: psi_lim, R_lim, Z_lim
integer :: i_elm_xpoint(2) 


character(len=MPI_MAX_PROCESSOR_NAME) :: name

logical, save :: axis_is_min

write(*,*) '***************************************'
write(*,*) '* JOREK2_diagno                       *'
write(*,*) '***************************************'

!my_id=0

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread

  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size

  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank

  ! --- Process command line arguments
  if ( my_id == 0 ) call jorek2help(n_cpu, nbthreads)

  CALL MPI_GET_PROCESSOR_NAME (name,resultlength,ierr)
  write(*,'(A,I5,2A)') '#MPI id, ProcessorName ', rank, ': ', name

  ! --- Initialize mode and mode_type arrays
  call det_modes()

  ! --- Remove file STOP_NOW if it exists
  if ( my_id == 0 ) then
    open(42, file='STOP_NOW', iostat=ierr)
    if ( ierr == 0 ) close(42, status='delete')
  end if


!call initialise_parameters(my_id, "__NO_FILENAME__")
call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")


! --- Write out all parameters defined in parameters and the namelist input file.
call log_parameters(my_id)

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis()                              ! define the basis functions at the Gaussian points

call broadcast_elements(my_id, element_list)                ! elements
call broadcast_nodes(my_id, node_list)                      ! nodes
call broadcast_phys(my_id)                                  ! physics parameters


!------------------------------------------------- find x-point(s)
xcase = 1
if (xpoint) then
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
  psi_bnd = psi_xpoint(1)
  if( (xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) ) then
    psi_bnd = psi_xpoint(2)
  endif
else
  psi_bnd = 0.d0
endif

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

if (my_id .eq. 0 ) then
  write(*,*) ' xcase,1st x-point:R,Z,psi: ',xcase, R_xpoint(1),Z_xpoint(1),psi_xpoint(1),psi_bnd
!   write(*,*) ' PSI_XPOINT : ',psi_xpoint,i_elm_xpoint
  write(*,*) ' PSI_AXIS : ',psi_axis,i_elm_axis
  write(*,*) ' RZ_AXIS : ', R_axis, Z_axis
endif


axis_is_min = axis_is_psi_minimum(node_list, element_list, R_axis, Z_axis, psi_axis) !Dump call to initialize


open(20,file='energies.txt')

write(20,'(A,25(A11,i3.3))') '      i      time',('          M',n_period*((in-1)/2),in=1,n_tor,2), &
                                                 ('          K',n_period*((in-1)/2),in=1,n_tor,2)

do i=2,index_start

 !Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
 !            / (xtime(i)-xtime(i-1))
 !Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
 !            / (xtime(i)-xtime(i-1))

 !write(*,'(i7,f12.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i),growth_mag,growth_kin

 write(20,'(i7,f12.3,200e14.6)') i,xtime(i),energies(1,1,i),(energies(in,1,i)+energies(in+1,1,i),in=2,n_tor,2), &
                                            energies(1,2,i),(energies(in,2,i)+energies(in+1,2,i),in=2,n_tor,2)

enddo
close(20)

if (use_pellet) then

  open(20,file="pellet.txt")

  write(20,'(A,25(A11,i3.3))') '      i,      time,    pellet_R,    pellet_Z,   pellet_psi,  particles,   ablation'

  do i=1,index_start
    write(20,'(i7,f12.3,200e14.6)') i,xtime(i),xtime_pellet_R(i),xtime_pellet_Z(i),xtime_pellet_psi(i),xtime_pellet_particles(i), &
                                    xtime_phys_ablation(i)
  enddo
  close(20)

endif

if (using_spi) then

  open(20,file="abl_history.dat")

  write(20,'(A11)') 'time', 'total_abl_rate', 'total_abl_number'

  do i=1,index_start
    spi_abl_rate_tot = 0.0
    spi_abl_tot = 0.0
    do i_spi = 1, n_spi
      spi_abl_rate_tot = spi_abl_rate_tot + xtime_spi_ablation_rate(i_spi,i)
      spi_abl_tot = spi_abl_tot + xtime_spi_ablation(i_spi,i)
    end do
    write(20,'(i7,f12.3,200e14.6)') i,xtime(i), spi_abl_rate_tot, spi_abl_tot
  enddo
  close(20)

endif

!call populate_element_rtree(node_list, element_list)

call Integrals_3D(my_id,node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

!if (use_pellet) then
!   pellet_volume = total_pellet_volume
!   call update_pellet(my_id,node_list,element_list)
!end if

!------------------lowshape3bis outside
!Rplot(1) = 3.0
!Rplot(2) = 3.676
!Zplot(1) = -2.066
!Zplot(2) = -1.9265

!------------------ lowshape3bis inside
!Rplot(1) = 2.3213
!Rplot(2) = 2.8511
!Zplot(1) = -1.9178
!Zplot(2) = -2.0339

!------------------ lowshape3bis midplane
!Rplot(1) = 1.88
!Rplot(2) = 4.88
!Zplot(1) = 0.09
!Zplot(2) = 0.09

!-----------------lowshape7,8 (midplane)
!Rplot(1) = 1.9
!Rplot(2) = 4.2
!Zplot(1) = 0.07
!Zplot(2) = 0.07

!call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

!Rplot(1) = 1.0
!Rplot(2) = 3.5
!Zplot(1) = Z_axis
!Zplot(2) = Z_axis 

!call plot_profiles(node_list,element_list,Rplot,Zplot)


!call export_helena(node_list,element_list)

!----------------------------------------- plot profiles
!call begplt('profiles.ps')

!call plot_velocity_profile(node_list,element_list, 3.d0, 0.d0, 3.d0, 2.d0)

!call finplt

end program jorek2_diagno

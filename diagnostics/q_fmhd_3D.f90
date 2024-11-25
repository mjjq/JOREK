!> Compute q-profile for a fMHD JOREK restart file from all n magnetic field
program q_fmhd_3D

use constants
use mod_parameters
use data_structure
use nodes_elements, only:bnd_node_list,bnd_elm_list
use equil_info
use phys_module
use mod_import_restart
use mod_log_params
use mod_interp
use elements_nodes_neighbours
use mod_find_rz_nearby
use mpi_mod
use mod_element_rtree, only: populate_element_rtree
use mod_q_profile_fmhd_3D

implicit none

integer   :: my_id
integer   :: required,provided,StatInfo, n_cpu, ierr
integer*4 :: rank, comm_size 

real*8 :: RaxisMin, RaxisMax, deltaRrange
real*8 :: ZaxisMin, ZaxisMax, deltaZrange
integer, parameter :: NaxisR = 40
integer, parameter :: NaxisZ = 40
real*8 :: Rguess1D(NaxisR), Zguess1D(NaxisZ)
real*8 :: Rguess3D(NaxisR, NaxisZ), Zguess3D(NaxisR, NaxisZ)
real*8 :: distMissMatch3D(NaxisR, NaxisZ), distMissMatch3D_myid(NaxisR, NaxisZ)
integer :: iR, iZ, iZ_myid, NaxisZ_myid, min_i, min_j, Niterations, iit
real*8  :: fact_ne, fact_T, fact_V

real(8), allocatable :: R_poinc(:), Z_poinc(:)
integer :: count_poinc, count_poinc_stpts

character(len=512) :: strtmp, iphi_string, nphi_string
integer :: nr_stpts, ntour_stpts, nphi_stpts, curr_stpts, i_stpts
real*8  :: rr_stpts, zz_stpts, phi_stpts
integer :: iphi_out, jphi_start, i_local
integer, parameter :: nphi_out = 32
real*8, parameter  :: phi_start = 0.d0
real*8  :: phi_poinc_record(nphi_out + 1) ! the last surface is the same as the first one

real*8  :: avg_count, avg_sum, avg_val

required = MPI_THREAD_FUNNELED
call MPI_Init_thread(required, provided, StatInfo)
call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_cpu = comm_size

! --- Determine ID of each MPI proc
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank

if ( my_id == 0 ) then
  write(*,*) '***************************************'
  write(*,*) '* Calculate q-profile...              *'
  write(*,*) '***************************************'
end if

call det_modes()
call initialise_parameters(my_id,  "__NO_FILENAME__")
call log_parameters(my_id)
if ( my_id == 0 ) call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)
call broadcast_phys(my_id)  
call broadcast_elements(my_id, element_list)                ! elements
call broadcast_nodes(my_id, node_list)                      ! nodes
call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)
call populate_element_rtree(node_list, element_list)
call broadcast_equil_state(my_id)

#define READ_STPTS

#ifdef READ_STPTS

! --- Same 'stpts' file as jorek2_poincare.f90 is used
! --- only the first 3 columns (nr, R_start, Z_start) will be used here.
! --- Read start points from file 'stpts'.
!
! Example for a stpts file:
!   Eleven field lines will be started, the first ten between (1.7,0.0) and (1.8,0.0), the eleventh at (1.85,0.2)
! 
!   +-------------------------------------------------------
!   |# n_lines
!   |  11
!   |# nr   R_start   Z_start    psi_start   n_turns   n_phis
!   |   1    1.700      0.000     0.000      100       1500
!   |  10    1.800      0.000     0.000      200       1500
!   |  11    1.850      0.200     0.000      800       1500
!   +-------------------------------------------------------
!

if(my_id == 0) then
  open(21, file='stpts', status='old', action='read', iostat=ierr)

  if ( ierr .ne. 0 ) then ! stpts file don't exists, exit.
    write(*,*) 'ERROR, stpts file is not found!!'
    stop
  endif

  read(21, '(a)') strtmp ! read comment line (ignored)
  read(21,*) count_poinc_stpts
  if ( count_poinc_stpts < 1 ) then
    write(*,*) 'ERROR in stpts file: count_poinc must be >= 1.'
    stop
  end if
  read(21, '(A)') strtmp ! read comment line (ignored)

  count_poinc = count_poinc_stpts + count_poinc_stpts / 10 ! add additional 10% lines near the axis

  allocate(R_poinc(count_poinc))
  allocate(Z_poinc(count_poinc))
    
  curr_stpts = 0
  do
    if ( curr_stpts >= count_poinc_stpts ) exit
      
    read(21, *) nr_stpts, rr_stpts, zz_stpts, phi_stpts, ntour_stpts, nphi_stpts
      
    if ( ( nr_stpts == 1 ) .and. ( curr_stpts == 0 ) ) then
      R_poinc(1) = rr_stpts
      Z_poinc(1) = zz_stpts
    else if ( curr_stpts == 0 ) then
      write(*,*) 'ERROR in stpts file: first start point must be nr_stpts=1.'
      stop
    else if ( ( nr_stpts < 1 ) .or. ( nr_stpts > count_poinc_stpts ) ) then
      write(*,*) 'ERROR in stpts file: nr_stpts must be > 0 and < count_poinc_stpts.'
      stop
    else if ( nr_stpts <= curr_stpts ) then
      write(*,*) 'ERROR in stpts file: start points must be sorted in ascending nr_stpts-order.'
      stop
    else
        
      do i_stpts = curr_stpts + 1, nr_stpts
        R_poinc(i_stpts) = R_poinc(curr_stpts) + ( rr_stpts - R_poinc(curr_stpts) ) * ( real(i_stpts-curr_stpts) / real(nr_stpts-curr_stpts) )
        Z_poinc(i_stpts) = Z_poinc(curr_stpts) + ( zz_stpts - Z_poinc(curr_stpts) ) * ( real(i_stpts-curr_stpts) / real(nr_stpts-curr_stpts) )
      end do
        
    end if
      
    curr_stpts = nr_stpts
      
  end do
    
  close(21)

  call linspace(5e-3, 1.e-1, count_poinc - count_poinc_stpts, R_poinc(count_poinc_stpts + 1:count_poinc))
  call linspace(5e-3, 1.e-1, count_poinc - count_poinc_stpts, Z_poinc(count_poinc_stpts + 1:count_poinc))

endif

call MPI_BCAST(count_poinc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(count_poinc_stpts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call construct_all(count_poinc, nphi_out)

if(my_id > 0) then
  allocate(R_poinc(count_poinc))
  allocate(Z_poinc(count_poinc))
endif

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_BCAST(R_poinc,count_poinc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(Z_poinc,count_poinc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! --- Define starting points for field lines along outer midplane
responsible = .false.
do i = 1, npoints
  if ( ( real(my_id)/real(n_cpu)*npoints < i ) .and. ( real(my_id+1)/real(n_cpu)*npoints >= i ) ) then
    responsible(i) = .true.
    Rstart(i) = R_poinc(i)
    Zstart(i) = Z_poinc(i)
  end if
end do

deallocate(R_poinc)
deallocate(Z_poinc)

#else

call construct_all(node_list%n_nodes - 200, nphi_out)
responsible = .false.
do i = 1, npoints
  if ( ( real(my_id)/real(n_cpu)*npoints < i ) .and. ( real(my_id+1)/real(n_cpu)*npoints >= i ) ) then
    responsible(i) = .true.
    Rstart(i) = node_list%node(i+100)%x(1,1,1)
    Zstart(i) = node_list%node(i+100)%x(1,1,2)
  end if
end do
#endif

call MPI_Barrier(MPI_COMM_WORLD,ierr)

! --- Write out distribution of field lines among tasks
do i = 0, n_cpu-1
  if ( my_id == i ) then
    write(*,*) 'Task ', my_id, 'responsible for (field line number & radius):'
    do j = 1, npoints
      if (responsible(j)) write(*,*) j, Rstart(j), Zstart(j)
    end do
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
end do

! --- specify the range for searching magnetic axis
! --- Example of 'axis_range.txt'
!  +-------------------------------------------------------
!  | ! RaxisMin  RaxisMax  ZaxisMin  ZaxisMax  Niterations
!  | 1.69      1.89      -0.07     0.2       1
!  +-------------------------------------------------------

! --- Find magnetic axis in 3D 1/1 mode
if ( my_id == 0 ) then
  open(911, file='axis_range.txt', status='old', action='read', form='formatted')
  read(911, '(a)') strtmp ! read comment line (ignored)
  read(911, *) RaxisMin, RaxisMax, ZaxisMin, ZaxisMax, Niterations
  close(911)
  write(*,*) "search magnetic axis within: " 
  write(*,*) "R: ", RaxisMin, RaxisMax
  write(*,*) "Z: ", ZaxisMin, ZaxisMax
  write(*,*) "Niterations: ", Niterations
endif
call MPI_BCAST(RaxisMin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(RaxisMax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ZaxisMin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ZaxisMax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(Niterations,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

do iit = 1,Niterations
  distMissMatch3D = 0.d0
  distMissMatch3D_myid = 0.d0
  call linspace(RaxisMin, RaxisMax, NaxisR, Rguess1D)
  call linspace(ZaxisMin, ZaxisMax, NaxisZ, Zguess1D)
  call meshgrid(Rguess1D, Zguess1D, NaxisR, NaxisZ, Rguess3D, Zguess3D)

  NaxisZ_myid = int(NaxisZ/n_cpu)
  do iZ_myid = 1, NaxisZ_myid
    iZ = my_id * NaxisZ_myid + iZ_myid
    if(iZ > NaxisZ) cycle
    if(my_id == 0) write(*,*) "... iZ: ", iZ_myid, ', of ', NaxisZ_myid, ", iit: ", iit, ", of ", Niterations
    do iR = 1, NaxisR
      call find_axis_3D(node_list, element_list, Rguess3D(iR, iZ), Zguess3D(iR, iZ), distMissMatch3D_myid(iR, iZ))
    enddo
  enddo
  call MPI_AllReduce(distMissMatch3D_myid, distMissMatch3D, NaxisZ * NaxisR, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call find_min_index(distMissMatch3D, NaxisR, NaxisZ, min_i, min_j)
  if(my_id == 0) then
    if(iit == 1) then
      call save_matrix_as_ascii(Rguess3D, NaxisR, NaxisZ, "Rguess3D.txt")
      call save_matrix_as_ascii(Zguess3D, NaxisR, NaxisZ, "Zguess3D.txt")
      call save_matrix_as_ascii(distMissMatch3D, NaxisR, NaxisZ, "distMissMatch3D.txt")
    endif
    write(*,*) "find axis at phi = 0 plane @ (", Rguess3D(min_i, min_j), ", ", Zguess3D(min_i, min_j), ")"
    write(*,*) "Error = ", distMissMatch3D(min_i, min_j)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  RaxisMin = Rguess3D(min_i, min_j) - (RaxisMax - RaxisMin) / 4.d0
  RaxisMax = Rguess3D(min_i, min_j) + (RaxisMax - RaxisMin) / 4.d0
  ZaxisMin = Zguess3D(min_i, min_j) - (ZaxisMax - ZaxisMin) / 4.d0
  ZaxisMax = Zguess3D(min_i, min_j) + (ZaxisMax - ZaxisMin) / 4.d0
enddo

if(my_id == 0) call record_axis_3D(node_list, element_list, Rguess3D(min_i, min_j), Zguess3D(min_i, min_j), stepPhiMax, Raxis3D, Zaxis3D)

call MPI_BCAST(Raxis3D,stepPhiMax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(Zaxis3D,stepPhiMax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

#ifdef READ_STPTS
jphi_start = nint(phi_start / stepsizePhi) + 1

do i_local = count_poinc_stpts + 1, npoints
  if(responsible(i_local)) then
    Rstart(i_local) = Rstart(i_local) + Raxis3D(jphi_start)
    Zstart(i_local) = Zstart(i_local) + Zaxis3D(jphi_start)
  endif
enddo
#endif

! --- Trace field lines
if ( my_id == 0 ) write(*,*) '*** ...start tracing... ***'  
call linspace(0.d0, 2.d0*PI, nphi_out + 1, phi_poinc_record)
call trance_field_lines_fmhd_3D(node_list, element_list, phi_start, phi_poinc_record(1:nphi_out), nphi_out)

fact_ne      = central_density * 1.d20 !* central_mass*MASS_PROTON    ! factor for rho
fact_T       = 1.d0 / ( MU_zero * central_density * 1.d20 * EL_CHG ) ! factor for T
fact_V       = 1.d0 / sqrt(MU_ZERO * central_density * 1.d20 * MASS_PROTON * central_mass)
ne_arr_poinc = ne_arr_poinc * fact_ne
Ti_arr_poinc = Ti_arr_poinc * fact_T
Te_arr_poinc = Te_arr_poinc * fact_T
VxBpara_poinc= VxBpara_poinc* fact_V

call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(my_id == 0) write(*,*) "Finished all tracing, start MPI_AllReduce..."

call MPI_AllReduce(MPI_IN_PLACE, psi_arr_poinc, num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, R_arr_poinc,   num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, Z_arr_poinc,   num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, q_arr_poinc,   num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, ne_arr_poinc,  num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, Te_arr_poinc,  num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, Ti_arr_poinc,  num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(MPI_IN_PLACE, VxBpara_poinc, num_poinc_max*npoints*nphi_out,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)


if(my_id == 0) then

  do i = 1, npoints
    do iphi_out = 1, nphi_out
      do j = 1, num_poinc_max
        if(R_arr_poinc(j,i,iphi_out) .le. 1.d-5 ) exit
      enddo
      num_poinc_min(i) = min(num_poinc_min(i), j - 1)
    enddo
    write(*,*) "Number of points for line ", i, ' is ', num_poinc_min(i)
  enddo

  write(nphi_string, '(I0)') nphi_out

  ! average the dynamo along the field line
  do i = 1, npoints
    avg_sum = 0.d0
    avg_count = 0.d0
    do iphi_out = 1, nphi_out
      do j = 1, num_poinc_min(i)
        avg_sum = avg_sum + VxBpara_poinc(j,i,iphi_out)
        avg_count = avg_count + 1.d0
      enddo
    enddo

    avg_val = avg_sum / avg_count
    do iphi_out = 1, nphi_out
      do j = 1, num_poinc_min(i)
        VxBpara_poinc(j,i,iphi_out) = avg_val
      enddo
    enddo
  enddo



  do iphi_out = 1, nphi_out

    write(*,*) 'Exporting No. ', iphi_out, ' plane of ', nphi_out
    write(iphi_string, '(I0)') iphi_out
    ! --- Write out result
    open(991, file='q_profile_3D_ip_'//trim(iphi_string)//'.dat', action='write', status='replace')
    write(991, *) '  Results @ phi = 2 * PI * ('//trim(iphi_string)//' - 1)/'//trim(nphi_string), ' = ', phi_poinc_record(iphi_out)
    write(991,'(A5, 8A23)') 'id', 'Psi_N', 'R', 'Z', 'q', 'density (m-3)', 'Ti (eV)', 'Te(eV)', 'VxBpara(V/m)'
    do i = 1, npoints
      do j = 1, num_poinc_min(i)
        write(991,'(I5, 8es23.8)') i, psi_arr_poinc(j,i,iphi_out), R_arr_poinc(j,i,iphi_out), Z_arr_poinc(j,i,iphi_out), q_arr_poinc(j,i,iphi_out), ne_arr_poinc(j,i,iphi_out), Ti_arr_poinc(j,i,iphi_out), Te_arr_poinc(j,i,iphi_out), VxBpara_poinc(j,i,iphi_out)
      end do
    end do
    close(991)
    write(*,*) '***************************************'
    write(*,*) '* Exported in q_profile_3D_ip_'//trim(iphi_string)//'.dat  *'
    write(*,*) '***************************************'
  enddo
endif


call deconstruct_all()

call MPI_FINALIZE(ierr)

end program q_fmhd_3D

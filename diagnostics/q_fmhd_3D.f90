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

real(8), allocatable :: R_poinc(:), Z_poinc(:)
integer :: count_poinc

character(len=512) :: strtmp
integer :: ios

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

!---> Read poincare points
count_poinc = 0
if(my_id == 0) then
  ! Call the subroutine to read the data
  call read_poinc_RZ('poinc_R-Z.dat', R_poinc, Z_poinc, count_poinc)
endif

call MPI_BCAST(count_poinc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!---> Distribute the points
if(count_poinc > 0) then
  call construct_all(count_poinc)

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

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  deallocate(R_poinc)
  deallocate(Z_poinc)

else

  call construct_all(node_list%n_nodes - 200)
  ! call construct_all(40)

  ! --- Define starting points for field lines along outer midplane
  responsible = .false.
  do i = 1, npoints
    
    if ( ( real(my_id)/real(n_cpu)*npoints < i ) .and. ( real(my_id+1)/real(n_cpu)*npoints >= i ) ) then
      responsible(i) = .true.
      Rstart(i) = node_list%node(i+100)%x(1,1,1)
      Zstart(i) = node_list%node(i+100)%x(1,1,2)
      ! Rstart(i) = 1.7686100
      ! Zstart(i) = 0.0470903
    end if
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

endif

! ! --- Write out distribution of field lines among tasks
! do i = 0, n_cpu-1
!   if ( my_id == i ) then
!     write(*,*) 'Task ', my_id, 'responsible for (field line number & radius):'
!     do j = 1, npoints
!       if (responsible(j)) write(*,*) j, Rstart(j), Zstart(j)
!     end do
!   end if
!   call MPI_Barrier(MPI_COMM_WORLD,ierr)
! end do

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
    if(my_id == 0) write(*,*) "iZ: ", iZ_myid, ', of ', NaxisZ_myid, ", iit: ", iit, ", of ", Niterations
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
    write(*,*) "find_axis_3D done! The axis at phi = 0 plane is "
    write(*,*) "(R, Z) = (", Rguess3D(min_i, min_j), ", ", Zguess3D(min_i, min_j), ")"
    write(*,*) "With Error = ", distMissMatch3D(min_i, min_j)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  RaxisMin = Rguess3D(min_i, min_j) - (RaxisMax - RaxisMin) / 4.d0
  RaxisMax = Rguess3D(min_i, min_j) + (RaxisMax - RaxisMin) / 4.d0
  ZaxisMin = Zguess3D(min_i, min_j) - (ZaxisMax - ZaxisMin) / 4.d0
  ZaxisMax = Zguess3D(min_i, min_j) + (ZaxisMax - ZaxisMin) / 4.d0
enddo

call record_axis_3D(node_list, element_list, Rguess3D(min_i, min_j), Zguess3D(min_i, min_j), stepPhiMax, Raxis3D, Zaxis3D, my_id)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if ( my_id == 0 ) write(*,*) '*** ...start tracing... ***'

! --- Trace field lines
call trance_field_lines_fmhd_3D(node_list, element_list)

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_Reduce(Rstart,    Rstart_all,    npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(Zstart,    Zstart_all,    npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(psin_arr,  psin_arr_tot,  npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(phi_arr,   phi_arr_tot,   npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(R_arr,     R_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(B_p_intg_arr,B_p_intg_arr_tot,npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(psi_arr,   psi_arr_tot,   npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(F_arr,     F_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

! --- Write out result
if ( my_id == 0 ) then
  open(99, file='q_profile_3D.dat', action='write', status='replace')
  write(99,*) '# Psi_N | R | Z | q | R at LFS midplane [m] | ∫dl/(R^2Bp)/2Pi | Psi | F'
  do i = 1, npoints
    write(99,'(8es23.5)') psin_arr_tot(i), Rstart_all(i), Zstart_all(i), phi_arr_tot(i)/(2.d0*PI), R_arr_tot(i), B_p_intg_arr_tot(i), psi_arr_tot(i), F_arr_tot(i)
  end do
  close(99)
  write(*,*) '***************************************'
  write(*,*) '* q-profile written in q_profile_3D.dat  *'
  write(*,*) '***************************************'
end if

call deconstruct_all()

call MPI_FINALIZE(ierr)

end program q_fmhd_3D

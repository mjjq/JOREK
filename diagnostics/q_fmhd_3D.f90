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
real*8  :: fact_ne, fact_T

real(8), allocatable :: R_poinc(:), Z_poinc(:)
integer :: count_poinc, count_poinc_stpts

character(len=512) :: strtmp
integer :: ios
integer :: nr_stpts, ntour_stpts, nphi_stpts, curr_stpts, i_stpts
real*8  :: rr_stpts, zz_stpts, phi_stpts
real*8  :: phi_start
integer :: iphi_out, nphi_out, jphi_start, i_local
character(len=100) :: iphi_string

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

#else

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
#endif

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

fact_ne      = central_density * 1.d20 !* central_mass*MASS_PROTON    ! factor for rho
fact_T       = 1.d0 / ( MU_zero * central_density * 1.d20 * EL_CHG ) ! factor for T

if ( my_id == 0 ) write(*,*) '*** ...start tracing... ***'

nphi_out = 32
do iphi_out = 1, nphi_out

  write(*,*) 'Tracing No. ', iphi_out, ' plane of ', nphi_out
  phi_start = 2 * PI / nphi_out * (iphi_out - 1)
  write(iphi_string, '(I0)') iphi_out

#ifdef READ_STPTS
  jphi_start = nint(phi_start / stepsizePhi) + 1

  do i_local = count_poinc_stpts + 1, npoints
    if(responsible(i_local)) then
      Rstart(i_local) = Rstart(i_local) + Raxis3D(jphi_start)
      Zstart(i_local) = Zstart(i_local) + Zaxis3D(jphi_start)
    endif
  enddo
#endif

  ! phi_start = 0.d0
  ! --- Trace field lines
  call trance_field_lines_fmhd_3D(node_list, element_list, phi_start)

  ! ne_arr = ne_arr * fact_ne
  ! Ti_arr = Ti_arr * fact_T
  ! Te_arr = Te_arr * fact_T

  ne_arr_poinc = ne_arr_poinc * fact_ne
  Ti_arr_poinc = Ti_arr_poinc * fact_T
  Te_arr_poinc = Te_arr_poinc * fact_T

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! call MPI_Reduce(Rstart,    Rstart_all,    npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(Zstart,    Zstart_all,    npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(psin_arr,  psin_arr_tot,  npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(phi_arr,   phi_arr_tot,   npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(R_arr,     R_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(B_p_intg_arr,B_p_intg_arr_tot,npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(psi_arr,   psi_arr_tot,   npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(F_arr,     F_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  ! call MPI_Reduce(ne_arr,   ne_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(Ti_arr,   Ti_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  ! call MPI_Reduce(Te_arr,   Te_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  call MPI_Reduce(psi_arr_poinc, psi_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(R_arr_poinc,     R_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(Z_arr_poinc,     Z_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(q_arr_poinc,     q_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(ne_arr_poinc,   ne_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(Te_arr_poinc,   Te_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(Ti_arr_poinc,   Ti_arr_poinc_tot, num_poinc_max*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! ! --- Write out result
  ! if ( my_id == 0 ) then
  !   open(99, file='q_profile_3D_simp_ip'//trim(iphi_string)//'.dat', action='write', status='replace')
  !   ! write(99,*) '# Psi_N | R | Z | q | R at LFS midplane [m] | ∫dl/(R^2Bp)/2Pi | Psi | F'
  !   write(99,*) '# id | Psi_N | R | Z | q | density (m-3) | Ti (eV) | Te(eV)'
  !   do i = 1, npoints
  !     ! write(99,'(8es23.5)') psin_arr_tot(i), Rstart_all(i), Zstart_all(i), phi_arr_tot(i)/(2.d0*PI), R_arr_tot(i), B_p_intg_arr_tot(i), psi_arr_tot(i), F_arr_tot(i)
  !     write(99,'(I5, 7es23.8)') i, psin_arr_tot(i), Rstart_all(i), Zstart_all(i), phi_arr_tot(i)/(2.d0*PI), ne_arr_tot(i), Ti_arr_tot(i), Te_arr_tot(i)
  !   end do
  !   close(99)
  !   write(*,*) '***************************************'
  !   write(*,*) '* q-profile written in q_profile_3D_simp.dat  *'
  !   write(*,*) '***************************************'
  ! end if

  ! --- Write out result
  if ( my_id == 0 ) then
    open(991, file='q_profile_3D_newip_'//trim(iphi_string)//'.dat', action='write', status='replace')
    write(991,*) '# id | Psi_N | R | Z | q | density (m-3) | Ti (eV) | Te(eV)'
    do i = 1, npoints
      do j = 1, num_poinc_max
          if (R_arr_poinc_tot(j,i).le.1.d-5) then 
            write(991,*) '    ' ! new line
            exit
          endif
          write(991,'(I5, 7es23.8)') i, psi_arr_poinc_tot(j,i), R_arr_poinc_tot(j,i), Z_arr_poinc_tot(j,i), q_arr_poinc_tot(j,i), ne_arr_poinc_tot(j,i), Ti_arr_poinc_tot(j,i), Te_arr_poinc_tot(j,i)
      end do
    end do
    close(991)
    write(*,*) '***************************************'
    write(*,*) '* q-profile written in q_profile_3D.dat  *'
    write(*,*) '***************************************'
  end if

#ifdef READ_STPTS
  jphi_start = nint(phi_start / stepsizePhi) + 1

  do i_local = count_poinc_stpts + 1, npoints
    if(responsible(i_local)) then
      Rstart(i_local) = Rstart(i_local) - Raxis3D(jphi_start)
      Zstart(i_local) = Zstart(i_local) - Zaxis3D(jphi_start)
    endif
  enddo
#endif

enddo

call deconstruct_all()

call MPI_FINALIZE(ierr)

end program q_fmhd_3D

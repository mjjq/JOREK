!> Compute q-profile for a fMHD JOREK restart file from n=0 magnetic field
program q_fmhd

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
use mod_q_profile_fmhd

implicit none

integer   :: my_id
integer   :: required,provided,StatInfo, n_cpu, ierr
integer*4 :: rank, comm_size 
real*8    :: psin_arr_tot(npoints) = 0.d0, phi_arr_tot(npoints) = 0.d0, R_arr_tot(npoints) = 0.d0
! --- Used to integrate 1/(2*pi)*int(dl/(R^2 B_poloidal))
real*8    :: B_p_intg_arr_tot(npoints) = 0.d0, psi_arr_tot(npoints) = 0.d0
real*8    :: F_arr_tot(npoints) = 0.d0

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

if ( my_id == 0 ) write(*,*) '*** ...start tracing... ***'

! --- Define starting points for field lines along outer midplane
responsible = .false.
do i = 1, npoints
  
  if ( ( real(my_id)/real(n_cpu)*npoints < i ) .and. ( real(my_id+1)/real(n_cpu)*npoints >= i ) ) then
    responsible(i) = .true.
    Rstart(i) = ES%R_axis + REAL(i) * (ES%R_midpl(2) - ES%R_axis - 0.002d0)/npoints
    Zstart(i) = ES%Z_axis
  end if
  
end do

! --- Write out distribution of field lines among tasks
do i = 0, n_cpu-1
  if ( my_id == i ) then
    write(*,*) 'Task ', my_id, 'responsible for (field line number & radius):'
    do j = 1, npoints
      if (responsible(j)) write(*,*) j, Rstart(j)
    end do
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
end do

! --- Trace field lines
call trance_field_lines_fmhd(node_list, element_list)

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_Reduce(psin_arr,  psin_arr_tot,  npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(phi_arr,   phi_arr_tot,   npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(R_arr,     R_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(B_p_intg_arr,B_p_intg_arr_tot,npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(psi_arr,   psi_arr_tot,   npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(F_arr,     F_arr_tot,     npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

! --- Write out result
if ( my_id == 0 ) then
  open(99, file='q_profile.dat', action='write', status='replace')
  write(99,*) '# Psi_N | Safety factor | R at LFS midplane [m] | ∫dl/(R^2Bp)/2Pi | Psi | F'
  do i = 1, npoints
    write(99,'(6es23.5)') psin_arr_tot(i), phi_arr_tot(i)/(2.d0*PI), R_arr_tot(i), B_p_intg_arr_tot(i), psi_arr_tot(i), F_arr_tot(i)
  end do
  close(99)
  write(*,*) '***************************************'
  write(*,*) '* q-profile written in q_profile.dat  *'
  write(*,*) '***************************************'
end if
call MPI_FINALIZE(ierr)

end program q_fmhd

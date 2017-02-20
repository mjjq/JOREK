!> Broadcast the numerical input profiles from MPI thread 0 to the others
subroutine broadcast_RMP_profiles(my_id, bnd_node_list)
use tr_module
use phys_module
use data_structure
use mpi_mod

implicit none

integer,                   intent(in) :: my_id
type (type_bnd_node_list), intent(in) :: bnd_node_list

integer :: ierr

if ( RMP_on ) then
  write(*,*) ' broadcast Number_RMP_harmonics=', Number_RMP_harmonics
  if ( my_id /= 0 ) then
    if (allocated(psi_RMP_cos))         call tr_deallocate(psi_RMP_cos,"psi_RMP_cos",CAT_UNKNOWN)
    if (allocated(dpsi_RMP_cos_dR))     call tr_deallocate(dpsi_RMP_cos_dR,"dpsi_RMP_cos_dR",CAT_UNKNOWN)
    if (allocated(dpsi_RMP_cos_dZ))     call tr_deallocate(dpsi_RMP_cos_dZ,"dpsi_RMP_cos_dZ",CAT_UNKNOWN)
    if (allocated(psi_RMP_sin))         call tr_deallocate(psi_RMP_sin,"psi_RMP_sin",CAT_UNKNOWN)
    if (allocated(dpsi_RMP_sin_dR))     call tr_deallocate(dpsi_RMP_sin_dR,"dpsi_RMP_sin_dR",CAT_UNKNOWN)
    if (allocated(dpsi_RMP_sin_dZ))     call tr_deallocate(dpsi_RMP_sin_dZ,"dpsi_RMP_sin_dZ",CAT_UNKNOWN)

    call tr_allocate(psi_RMP_cos,1,    bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_cos",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_cos_dR,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dR",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_cos_dZ,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dZ",CAT_UNKNOWN)
    call tr_allocate(psi_RMP_sin,1,    bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_sin",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_sin_dR,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dR",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_sin_dZ,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dZ",CAT_UNKNOWN)
  end if
  call MPI_BCAST(psi_RMP_cos,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dpsi_RMP_cos_dR,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dpsi_RMP_cos_dZ,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(psi_RMP_sin,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dpsi_RMP_sin_dR,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dpsi_RMP_sin_dZ,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
end if

return
end subroutine broadcast_RMP_profiles

module centralization_mod



implicit none



contains



!> Collect a distributed matrix on a single MPI process (the master task of the MPI group)
subroutine centralization_harmonic(my_id, my_id_n, n_cpu_n, MPI_COMM_N)
  use data_structure 
  use global_distributed_matrix
  use mumps_module  
  use mpi_mod
  
  implicit none
  
  ! --- Routine parameters
  integer, intent(in)          :: my_id, my_id_n, n_cpu_n, MPI_COMM_N
  
  ! --- Local variables
  integer, allocatable         :: nz_array(:), disp_array(:)
  integer                      :: nz_total, i, ierr

  call MPI_Barrier(MPI_COMM_WORLD,ierr) !###

    write(*,*) my_id, my_id_n, '###A0', n_cpu_n
    write(*,*) my_id, my_id_n, '###A0', allocated(nz_array), allocated(disp_array)
    call system ('sleep 5s')
    
  allocate(nz_array  (n_cpu_n))
    write(*,*) my_id, my_id_n, '###A0b'
    call system ('sleep 5s')
    
  allocate(disp_array(n_cpu_n))
  
    write(*,*) my_id, my_id_n, '###A1'
    call system ('sleep 5s')
    
  ! --- Determine matrix dimension n and number of nonzeros nz
  if (n_cpu_n > 1) then
    write(*,*) my_id, my_id_n, '###A2'
    call system ('sleep 5s')
    
    call MPI_GATHER(nz_glob_harm, 1, MPI_INTEGER, nz_array, 1, MPI_INTEGER, 0, MPI_COMM_N, ierr) 
    disp_array = 0 
    nz_total   = 0
    if (my_id_n .eq. 0) then
      do i = 2, n_cpu_n  
         disp_array(i) = disp_array(i-1) + nz_array(i-1)
      enddo  
      nz_total = disp_array(n_cpu_n) + nz_array(n_cpu_n) 
    endif 
    mumps_par%nz  = nz_total  !nz_glob_harm 
    mumps_par%n   = ndof_glob_harm 
  else
    write(*,*) my_id, my_id_n, '###A3'
    call system ('sleep 5s')
    
    mumps_par%nz  = nz_glob_harm
    mumps_par%n   = ndof_glob_harm
  endif 
    
    write(*,*) '###dim_nz ', my_id, my_id_n, mumps_par%nz
    write(*,*) '###dim_n  ', my_id, my_id_n, mumps_par%n
    write(*,*) '###dim_agh', my_id, my_id_n, size(A_glob_harm,1)
    write(*,*) '###dim_igh', my_id, my_id_n, size(irn_glob_harm,1)
    write(*,*) '###dim_jgh', my_id, my_id_n, size(jcn_glob_harm,1)
    write(*,*) '###dim_rgh', my_id, my_id_n, size(rhs_glob_harm,1)
    
    write(*,*) my_id, my_id_n, '###A'
    call system ('sleep 5s')
    
  ! --- Allocate arrays for centralized matrix
  if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,  "dh_mumps_par%A",  CAT_DMATRIX)
  if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
  if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)
  if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs",CAT_DMATRIX)

    write(*,*) my_id, my_id_n, '###B', mumps_par%n, mumps_par%nz
    call system ('sleep 5s')
    
  call tr_allocatep(mumps_par%A,  1,mumps_par%nz,"dh_mumps_par%A",  CAT_DMATRIX)
  call tr_allocatep(mumps_par%irn,1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
  call tr_allocatep(mumps_par%jcn,1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,1,mumps_par%n, "dh_mumps_par%rhs",CAT_DMATRIX)
  mumps_par%a = 42.d0 !###
  mumps_par%irn = 43  !###
  mumps_par%jcn = 44  !###

    write(*,*) my_id, my_id_n, '###C'
    call system ('sleep 5s')
    
  ! --- Centralize matrix (if it was not distributed, copy it into the right data structure)
  if (n_cpu_n > 1) then
    write(*,*) my_id, my_id_n, '###D'
    call system ('sleep 5s')
    call MPI_GATHERV(A_glob_harm, nz_glob_harm, MPI_DOUBLE_PRECISION, mumps_par%A, nz_array, disp_array,&
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_N, ierr)
    call MPI_GATHERV(irn_glob_harm, nz_glob_harm, MPI_INTEGER, mumps_par%irn, nz_array, disp_array,&
                     MPI_INTEGER, 0, MPI_COMM_N, ierr)
    call MPI_GATHERV(jcn_glob_harm, nz_glob_harm, MPI_INTEGER, mumps_par%jcn, nz_array, disp_array,&
                     MPI_INTEGER, 0, MPI_COMM_N, ierr)
  else 
    write(*,*) my_id, my_id_n, '###E1', size(mumps_par%irn,1), size(irn_glob_harm,1)
    call system ('sleep 5s')
    mumps_par%irn(1:mumps_par%nz) = irn_glob_harm(1:mumps_par%nz)
    write(*,*) my_id, my_id_n, '###E2', size(mumps_par%jcn,1), size(jcn_glob_harm,1)
    call system ('sleep 5s')
    mumps_par%jcn(1:mumps_par%nz) = jcn_glob_harm(1:mumps_par%nz)
    write(*,*) my_id, my_id_n, '###E', size(mumps_par%A,1), size(A_glob_harm,1)
    call system ('sleep 5s')
    mumps_par%A(1:mumps_par%nz)   = A_glob_harm(1:mumps_par%nz)
  endif
  
    write(*,*) my_id, my_id_n, '###F'
    call system ('sleep 5s')
  mumps_par%rhs = rhs_glob_harm

    write(*,*) my_id, my_id_n, '###G'
    call system ('sleep 5s')
  if ( allocated(nz_array) )   deallocate(nz_array)
  if ( allocated(disp_array) ) deallocate(disp_array)

    write(*,*) my_id, my_id_n, '###H'
    call system ('sleep 5s')
end subroutine centralization_harmonic

end module centralization_mod

module centralization_mod
  implicit none
contains  
subroutine centralization_harmonic(my_id_n, n_cpu_n, MPI_COMM_N)
  use data_structure 
  use global_distributed_matrix
  use mumps_module  
  use mpi_mod
  implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_bnd_element_list) :: bnd_elm_list
type (type_surface_list)     :: flux_list
type (type_element)          :: element
type (type_node)             :: nodes(n_vertex_max)



integer, intent(in)          :: my_id_n, n_cpu_n, MPI_COMM_N
integer, allocatable         :: nz_array(:)  
integer, allocatable         :: disp_array(:)  
integer                      :: nz_total, i, ierr     

      !----- collecting distributed matrix on a single MPI process
      if (.not.allocated(nz_array))   allocate(nz_array(n_cpu_n)) 
      if (.not.allocated(disp_array)) allocate(disp_array(n_cpu_n))
 
      call MPI_GATHER(nz_glob_harm, 1, MPI_INTEGER, nz_array, 1, MPI_INTEGER, 0, MPI_COMM_N, ierr) 

      disp_array = 0 
      nz_total   = 0
      if(my_id_n .eq. 0) then
        do i = 2, n_cpu_n  
           disp_array(i) = disp_array(i-1) + nz_array(i-1)
        enddo  
      nz_total = disp_array(n_cpu_n) + nz_array(n_cpu_n) 
      endif 
 
      mumps_par%nz = nz_total  
      mumps_par%n  = ndof_glob_harm

      if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"dh_mumps_par%A",CAT_DMATRIX)
      if (associated(mumps_par%irn))  call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
      if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)

      call tr_allocatep(mumps_par%A,1,mumps_par%nz,"dh_mumps_par%A",CAT_DMATRIX)
      call tr_allocatep(mumps_par%irn,1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
      call tr_allocatep(mumps_par%jcn,1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)

      if (associated(mumps_par%rhs))  call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs",CAT_DMATRIX)
      call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"dh_mumps_par%rhs",CAT_DMATRIX)

      call MPI_GATHERV(A_glob_harm, nz_glob_harm, MPI_DOUBLE_PRECISION, mumps_par%A, nz_array, disp_array,&
                       MPI_DOUBLE_PRECISION, 0, MPI_COMM_N, ierr)
      call MPI_GATHERV(irn_glob_harm, nz_glob_harm, MPI_INTEGER, mumps_par%irn, nz_array, disp_array,&
                       MPI_INTEGER, 0, MPI_COMM_N, ierr)
      call MPI_GATHERV(jcn_glob_harm, nz_glob_harm, MPI_INTEGER, mumps_par%jcn, nz_array, disp_array,&
                       MPI_INTEGER, 0, MPI_COMM_N, ierr)
      
      mumps_par%rhs = rhs_glob_harm

      if ( allocated(nz_array) )   deallocate(nz_array) 
      if ( allocated(disp_array) ) deallocate(disp_array) 

end subroutine centralization_harmonic
end module centralization_mod

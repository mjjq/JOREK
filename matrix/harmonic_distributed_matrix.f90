!> Contains data structures and some routines related to the global matrix and right hand sides.
!!
!! The sparse matrix and the corresponding right hand side represent the system of equations
!! given to the solver for calculating the time evolution of the physical quantities.
module harmonic_distributed_matrix
  
  implicit none
  
  public

  !use mumps_module    
  !real*8,  allocatable, target  :: mumps_par%A(:)    !< Distributed global matrix
  !real*8,  allocatable, target  :: mumps_par%rhs(:)  !< Distributed global right hand side
  !integer, allocatable, target  :: mumps_par%irn(:)  !< Row indices for coordinate format sparse matrix (or CSR)
  !integer, allocatable, target  :: mumps_par%jcn(:)  !< Column indices for coordinate format sparse matrix (or CSR)
  !real*8,  allocatable  :: mumps_par%A(:)    !< Distributed global matrix
  !real*8,  allocatable  :: mumps_par%rhs(:)  !< Distributed global right hand side
  !integer, allocatable  :: mumps_par%irn(:)  !< Row indices for coordinate format sparse matrix (or CSR)
  !integer, allocatable  :: mumps_par%jcn(:)  !< Column indices for coordinate format sparse matrix (or CSR)
  integer, allocatable :: ijA_index_harm(:,:), ijA_size_harm(:), irn_jcn_harm(:,:) !< contains the structure of the sparse matrix (to fill in CSR format)
  integer, allocatable  :: irn_glob_harm(:)  !< Row indices for coordinate format sparse matrix (or CSR)
  integer, allocatable  :: jcn_glob_harm(:)  !< Column indices for coordinate format sparse matrix (or CSR)

!  real*8,  allocatable :: deltas(:)                                 !< solution from previous step
!  real*8,  allocatable :: column_scaling(:)                         !< column scaling of the global matrix
 ! integer, allocatable :: local_index_start(:), local_index_end(:)  !< range of indices local to one MPI process 
!  integer              :: ndof_glob
  integer              :: n_matrix_block_size_harm                       !< Size of a matrix block (n_var x n_tor)
  integer              :: n_glob_harm, nz_glob_harm                       !< Size of a matrix block (n_var x n_tor)
  
  
  
end module harmonic_distributed_matrix

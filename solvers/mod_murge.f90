
module murge_module             ! Murge specific variables
  implicit none
  
  include "murge.inc"
  !include "hips.inc"
  ! Indicate which solver is used
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_solver  

  ! Solver identification number
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_id      
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_id_prod
  INTEGER                                     :: murge_harmonic
  ! Local number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_local_n 
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_local_n_prod 
  ! Global number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_global_n 
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_global_n_prod 
  ! Number of dof by node
  INTEGER(KIND=MURGE_INTS_KIND)               :: murge_ndof

  ! Local element list
  INTEGER(KIND=MURGE_INTS_KIND), allocatable  :: murge_loc2glob(:)
  INTEGER(KIND=MURGE_INTS_KIND), allocatable  :: murge_loc2glob_prod(:)
  ! Global element list
  INTEGER(KIND=MURGE_INTS_KIND), allocatable  :: murge_glob2loc(:)
  INTEGER(KIND=MURGE_INTS_KIND), allocatable  :: murge_glob2loc_prod(:)

  ! Indicate if we want to use murge or classical interface
  logical                                     :: use_murge

  ! Indicate if we want to use murge element building
  logical                                     :: use_murge_element
  logical                                     :: use_hips
  ! Indicate if murge has been initialized
  logical                                     :: murge_initialised

  ! Murge right-hand-side member
  integer(KIND=MURGE_INTS_KIND)               :: murge_rhs
  parameter (murge_rhs=0)

  ! Indicate if the matrix is symmetric
  integer(KIND=MURGE_INTS_KIND)               :: murge_sym
  parameter (murge_sym=MURGE_BOOLEAN_FALSE)

  ! Number of threads
  integer(KIND=MURGE_INTS_KIND)               :: murge_nthrd =1

  ! Number of iteration in refinement
  integer(KIND=MURGE_INTS_KIND)               :: murge_iter
  parameter (murge_iter=10)

  integer(KIND=MURGE_INTS_KIND)               :: murge_ricar
  parameter (murge_ricar=0)

  integer(KIND=MURGE_INTS_KIND)               :: murge_iluk
  parameter (murge_iluk=3)

  integer(KIND=MURGE_INTS_KIND)               :: murge_amalg
  parameter (murge_amalg=5)

  real*8                                      :: murge_epsilon
  parameter (murge_epsilon=1.d-12)

  real*8                                      :: murge_pivot
  parameter (murge_pivot=1.d-64)
  integer                                     :: murge_comm

  contains
    !
    ! Subroutine: murge_add_one_entry
    ! 
    ! Add one entry to the product and/or harminic matrix.
    ! 
    ! Parameters :
    !   index_node  - row node index
    !   k           - row var index
    !   in          - row tor index
    !   index_node2 - col node index
    !   k2          - col var index
    !   in2         - col tor index
    !   zbig        - value
    !   solve_only  - Do not add to harmonic matrix if .true.
    !   gmres       - Do not add to product matrix if .false.
    !
    subroutine murge_add_one_entry( index_node, k, in, index_node2, k2,&
         &                          in2, zbig, murge_ntor, solve_only, gmres )
      use parameters
      INTEGER :: index_node,  k,  in
      INTEGER :: index_node2, k2, in2
      REAL*8  :: zbig
      integer :: murge_ntor
      integer :: ierr
      logical :: solve_only, gmres
      integer(KIND=MURGE_INTS_KIND) :: row_idx, col_idx  

      if (.not. solve_only) then
         row_idx = murge_ndof * (index_node - 1) + (k -1)*murge_ntor + 1 
         col_idx = murge_ndof * (index_node2- 1) + (k2-1)*murge_ntor + 1
         if (in /= 1) row_idx = row_idx + mod(in, 2) 
         if (in2 /= 1) col_idx = col_idx + mod(in2, 2) 
         CALL MURGE_ASSEMBLYSETVALUE( murge_id, row_idx, col_idx, zbig,&
              &                       ierr )
         IF (ierr /= MURGE_SUCCESS) THEN
            WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE ", &
                 "I", index_node, murge_ndof, k, in, &
                 "J", index_node2, k2, in2
            STOP
         END IF
      end if
      if (gmres) then
         row_idx = n_tor*n_var * (index_node - 1) + (k -1)*n_tor + in
         col_idx = n_tor*n_var * (index_node2- 1) + (k2-1)*n_tor + in2 
         CALL MURGE_ASSEMBLYSETVALUE( murge_id_prod, row_idx, col_idx, zbig,&
              &                       ierr )
         IF (ierr /= MURGE_SUCCESS) THEN
            WRITE (*,*)  "ERROR in MURGE_ASSEMBLYSETVALUE(prod) ", &
                 "I", index_node, n_tor, n_var, k, in, &
                 "J", index_node2, k2, in2
            STOP
         END IF
      end if
    end subroutine murge_add_one_entry

endmodule murge_module

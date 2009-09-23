module murge_module             ! Murge specific variables
  
  include "murge.inc"
  ! Indicate which solver is used
  INTEGER(KIND=MURGE_INTS_KIND)               :: solver  
  parameter( solver = MURGE_SOLVER_PASTIX ) 

  ! Solver identification number
  INTEGER(KIND=MURGE_INTS_KIND)               :: id      

  ! Local number of element
  INTEGER(KIND=MURGE_INTS_KIND)               :: local_n 

  ! Local element list
  INTEGER(KIND=MURGE_INTS_KIND), allocatable  :: loc2glob(:)

  ! Indicate if we want to use murge or classical interface
  logical                                     :: use_murge

  ! Indicate if we want to use murge element building
  logical                                     :: use_murge_element

  ! Indicate if murge has been initialized
  logical                                     :: murge_initialised

  ! Murge right-hand-side member
  integer(KIND=MURGE_INTS_KIND)               :: murge_rhs
  parameter (murge_rhs=0)

  ! Indicate if the matrix is symmetric
  integer(KIND=MURGE_INTS_KIND)               :: murge_sym
  parameter (murge_sym=MURGE_BOOLEAN_FALSE)

  ! Number of threads
  integer(KIND=MURGE_INTS_KIND)               :: murge_nthrd
  parameter (murge_nthrd=1)

  ! Number of iteration in refinement
  integer(KIND=MURGE_INTS_KIND)               :: murge_iter
  parameter (murge_iter=2)

  integer(KIND=MURGE_INTS_KIND)               :: murge_ricar
  parameter (murge_ricar=0)

  integer(KIND=MURGE_INTS_KIND)               :: murge_iluk
  parameter (murge_iluk=3)

  integer(KIND=MURGE_INTS_KIND)               :: murge_amalg
  parameter (murge_amalg=5)

  real*8                                      :: murge_epsilon
  parameter (murge_epsilon=1.d-12)

  real*8                                      :: murge_pivot
  parameter (murge_pivot=1.d-32)

endmodule murge_module

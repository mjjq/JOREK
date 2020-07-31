!> Module enabling JOREK to use the Watson Sparse Matrix Package (WSMP) as a solver
!! for the time-evolution. Note, that it is currently not used for the equilibrium.
module wsmp_module
  
  implicit none
  
  private
  public PWGSMP__allocate, PWGSMP__deallocate, PWGSMP__initialize_matrix,                          &
    PWGSMP__initialize_solver, PWGSMP__LU_factorization, PWGSMP__back_substitution
  
  ! --- Hard-coded constants
  logical, parameter     :: PWSMP__copy_matrix = .false. !< Copy matrix instead of using pointers?
  logical, parameter     :: PWSMP__verbose     = .false. !< Output additional debug information?
  
  !> Data structure that contains the system of equations to solve and the WSMP parameters.
  type :: PWGSMP_STRUCT
    ! --- THE MATRIX (allocation/contents are valid on MPI rank 0 only) ---
    integer              :: rank      !< rank of the matrix (aka "N")
    integer              :: nnz       !< number of non-zero entries
    integer, pointer     :: ir(:)     !< row indices
    integer, pointer     :: jc(:)     !< column indices
    real*8,  pointer     :: m(:)      !< THE MATRIX itself
    real*8,  pointer     :: b(:)      !< the right hand side of size (ldb x nrhs), may be 0 iff N_i==0
    !--- auxiliary fields as documented by the manual for PWSSMP (p.41ff), must be initialized separately ---
    integer              :: n         !< number of rows of the RHS residing on process i ("N_i")
    integer              :: ldb       !< must be >= N_i
    integer              :: nrhs      !< must be the same on all proc's, typically "1"
    real*8, allocatable  :: rmsic(:)  !< must be of size (ldb x nrhs), may be 0 iff N_i==0
    integer              :: iparm(64) !< integer parameters controlling WSMP
    real*8               :: dparm(64) !< integer parameters controlling WSMP
  end type PWGSMP_STRUCT
  
  type(PWGSMP_STRUCT)    :: PWGSMP__matrix  !< WSMP data structure
  integer                :: ierr            !< Error code
  
  
  
  contains
  
  
  
  !> Allocate memory for the WSMP data structure
  subroutine PWGSMP__allocate(n, nnz, my_id_n)
    
    integer, intent(in) :: n, nnz, my_id_n
    
    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__allocate() ..."
    
    if (PWSMP__copy_matrix) then
      if (associated(PWGSMP__matrix%ir)) deallocate(PWGSMP__matrix%ir)
      if (associated(PWGSMP__matrix%jc)) deallocate(PWGSMP__matrix%jc)
      if (associated(PWGSMP__matrix%m))  deallocate(PWGSMP__matrix%m)
    endif
    
    if ( PWSMP__copy_matrix .and. (my_id_n .eq. 0) ) then
      allocate(PWGSMP__matrix%ir(n+1))
      allocate(PWGSMP__matrix%jc(nnz))
      allocate(PWGSMP__matrix%m(nnz))
      PWGSMP__matrix%ir(:) = -1
      PWGSMP__matrix%jc(:) = -1
      PWGSMP__matrix%m(:)  = 0.0
    else
      PWGSMP__matrix%ir => null()
      PWGSMP__matrix%jc => null()
      PWGSMP__matrix%m  => null()
    endif
    
    PWGSMP__matrix%rank    = -1
    PWGSMP__matrix%nnz     = -1
    
    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__allocate() ..."
    
  end subroutine PWGSMP__allocate
  
  
  
  !> Deallocate memory of the WSMP data structure
  subroutine PWGSMP__deallocate()
    
    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__deallocate() ..."
    
    if (PWSMP__copy_matrix) then
      if (associated(PWGSMP__matrix%ir))   deallocate(PWGSMP__matrix%ir)
      if (associated(PWGSMP__matrix%jc))   deallocate(PWGSMP__matrix%jc)
      if (associated(PWGSMP__matrix%m))    deallocate(PWGSMP__matrix%m)
      if (associated(PWGSMP__matrix%b))    deallocate(PWGSMP__matrix%b)
    endif
    if (allocated(PWGSMP__matrix%rmsic)) deallocate(PWGSMP__matrix%rmsic)
    
    PWGSMP__matrix%ir => null()
    PWGSMP__matrix%jc => null()
    PWGSMP__matrix%m  => null()
    PWGSMP__matrix%b  => null()
    
    PWGSMP__matrix%rank  = -1
    PWGSMP__matrix%nnz   = -1
    
    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__deallocate() ..."
    
  end subroutine PWGSMP__deallocate
  
  
  
  !> Initialize the matrix in the WSMP data structure by either copying the data or by
  !! using pointers to the original data.
  subroutine PWGSMP__initialize_matrix(n, nnz, a, ir, jc, my_id_n)
    
    integer, intent(in)  ::     n ! dimension
    integer, intent(in)  ::   nnz ! number of non-zeros
    ! input matrix
    real*8, pointer      ::  a(:) ! matrix elements, length: nnz
    integer, pointer     :: ir(:) ! row indices,     length: nnz
    integer, pointer     :: jc(:) ! column indices,  length: nnz
    integer, intent(in)  :: my_id_n

    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__initialize_matrix() ..."

    if (my_id_n .eq. 0) then
      if ( PWSMP__copy_matrix .and. &
          .not.(      (associated(PWGSMP__matrix%ir))   &
                 .and.(associated(PWGSMP__matrix%jc))   &
                 .and.(associated(PWGSMP__matrix%m )) ) &
         ) then
        call wsmp_error_handler('PWGSMP__prepare() called on unallocated memory')
      endif

      PWGSMP__matrix%rank      = n
      PWGSMP__matrix%nnz       = nnz

      if (PWSMP__copy_matrix) then
        PWGSMP__matrix%ir(1:n+1) = ir(1:n+1)
        PWGSMP__matrix%jc(1:nnz) = jc(1:nnz)
        PWGSMP__matrix%m(1:nnz)  = a(1:nnz)
      else
        ! associate pointers
        PWGSMP__matrix%ir => ir
        PWGSMP__matrix%jc => jc
        PWGSMP__matrix%m  => a
      endif
    endif

    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__initialize_matrix() ..."
    
  end subroutine PWGSMP__initialize_matrix
  
  
  
  !> Initialize the WSMP solver.
  subroutine PWGSMP__initialize_solver(my_id_n, MPI_COMM_N)

    use mpi_mod
    !$ use omp_lib

    integer, intent(in) :: my_id_n
    integer, intent(in) :: MPI_COMM_N

    integer             :: PWGSMP_nthrd

    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__initialize_solver() ..."

    ! --- initialize task and thread parallelism
#ifdef USE_WSMP
    call wsetmpicomm(MPI_COMM_N)
#endif
#ifdef _OPENMP
    !$omp parallel default(none) shared(PWGSMP_nthrd)
    PWGSMP_nthrd = omp_get_num_threads()
    !$omp end parallel
#else
    PWGSMP_nthrd = 1
#endif
#ifdef USE_WSMP
    call wsetmaxthrds(PWGSMP_nthrd)
#endif


    call MPI_BCAST(PWGSMP__matrix%rank, 1, MPI_INTEGER, 0, MPI_COMM_N, ierr)

    ! --- prepare all variables to run PWGSMP in 0-master mode
    !     deallocate any allocated fields
    if (associated(PWGSMP__matrix%b))    deallocate(PWGSMP__matrix%b)
    if (allocated(PWGSMP__matrix%rmsic)) deallocate(PWGSMP__matrix%rmsic)
    
    ! --- allocate/fill fields which differ between MPI ranks
    if (my_id_n == 0) then
      PWGSMP__matrix%ldb = PWGSMP__matrix%n
      PWGSMP__matrix%n = PWGSMP__matrix%rank
      allocate(PWGSMP__matrix%rmsic(PWGSMP__matrix%n))
      allocate(PWGSMP__matrix%b(PWGSMP__matrix%n))
    else
      if (associated(PWGSMP__matrix%ir)) deallocate(PWGSMP__matrix%ir)
      if (associated(PWGSMP__matrix%jc)) deallocate(PWGSMP__matrix%jc)
      if (associated(PWGSMP__matrix%m))  deallocate(PWGSMP__matrix%m)
      PWGSMP__matrix%ir => null()
      PWGSMP__matrix%jc => null()
      PWGSMP__matrix%m  => null()
      
      PWGSMP__matrix%n = 0
      allocate(PWGSMP__matrix%rmsic(1))
      PWGSMP__matrix%b => null()
      PWGSMP__matrix%ldb = 1
    endif
    
    ! --- allocate/fill fields which are the same on all MPI ranks
    PWGSMP__matrix%rmsic(:) = 0.0
    PWGSMP__matrix%nrhs = 1
    PWGSMP__matrix%iparm(:) = 0
    PWGSMP__matrix%dparm(:) = 0.0

    ! --- fill iparm and dparm with default values
    PWGSMP__matrix%iparm(1) = 0
    PWGSMP__matrix%iparm(2) = 0
    PWGSMP__matrix%iparm(3) = 0

    PWGSMP__matrix%iparm(4)=1

    if (my_id_n .eq. 0) then
#ifdef USE_WSMP
      call ws_sortindices_d( PWGSMP__matrix%n,  PWGSMP__matrix%n,&
                             PWGSMP__matrix%ir, PWGSMP__matrix%jc,&
                             PWGSMP__matrix%m,  ierr )
#endif
      if (ierr .ne. 0) call wsmp_error_handler("PWGSMP index sorting failed.")
    endif

    if (my_id_n .eq. 0) then
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     PWGSMP__matrix%ir,    PWGSMP__matrix%jc,   PWGSMP__matrix%m,&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, null(),&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    else
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     null(),               null(),              null(),&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, null(),&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    endif
    if (PWGSMP__matrix%iparm(64) .ne. 0) call wsmp_error_handler("PWGSMP initialization failed.")

    PWGSMP__matrix%iparm(4)=1

    ! --- analysis ---
    PWGSMP__matrix%iparm(2) = 1
    PWGSMP__matrix%iparm(3) = 1

    if (my_id_n .eq. 0) then
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     PWGSMP__matrix%ir,    PWGSMP__matrix%jc,   PWGSMP__matrix%m,&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, null(),&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    else
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     null(),               null(),              null(),&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, null(),&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    endif

    if (PWGSMP__matrix%iparm(64) .ne. 0) then
      write(*,*) "PWGSMP__matrix%iparm(64): ", PWGSMP__matrix%iparm(64)
      call wsmp_error_handler("PWGSMP analysis failed.")
    endif

    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__initialize_solver() ..."

  end subroutine PWGSMP__initialize_solver
  
  
  
  !> Perform the LU factorization.
  subroutine PWGSMP__LU_factorization(my_id_n)
    use mpi_mod
    integer, intent(in)  :: my_id_n


    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__LU_factorization() ..."

    PWGSMP__matrix%iparm(2) = 2
    PWGSMP__matrix%iparm(3) = 2

    if (my_id_n .eq. 0) then
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     PWGSMP__matrix%ir,    PWGSMP__matrix%jc,   PWGSMP__matrix%m,&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, PWGSMP__matrix%rmsic,&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    else
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     null(),               null(),              null(),&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, PWGSMP__matrix%rmsic,&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    endif

    if (PWGSMP__matrix%iparm(64) .ne. 0) call wsmp_error_handler("PWGSMP symbolic factorization failed.")
    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__LU_factorization() ..."

  end subroutine PWGSMP__LU_factorization
  
  
  
  !> Perform the back substitution.
  subroutine PWGSMP__back_substitution(rhs, my_id_n)
    use mpi_mod
    real*8, pointer     :: rhs(:)
    integer, intent(in) :: my_id_n
    integer             :: i


    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__back_substitution() ..."

    if (my_id_n .eq. 0) then
      if (PWSMP__copy_matrix) then
        PWGSMP__matrix%b(1:PWGSMP__matrix%n) = rhs(1:PWGSMP__matrix%n)
      else
        PWGSMP__matrix%b => rhs
      endif
    endif

    PWGSMP__matrix%iparm(2) = 3
    PWGSMP__matrix%iparm(3) = 3

    if (my_id_n .eq. 0) then
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     PWGSMP__matrix%ir,    PWGSMP__matrix%jc,   PWGSMP__matrix%m,&
                  PWGSMP__matrix%b,     PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, PWGSMP__matrix%rmsic,&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    else
#ifdef USE_WSMP
      call pwgsmp(PWGSMP__matrix%n,     null(),               null(),              null(),&
                  null(),               PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, PWGSMP__matrix%rmsic,&
                  PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
#endif
    endif

    if (PWGSMP__matrix%iparm(64) .ne. 0) call wsmp_error_handler("PWGSMP back substitution failed.")

    if (my_id_n .eq. 0) then
      if (PWSMP__copy_matrix) then
        rhs(1:PWGSMP__matrix%n) = PWGSMP__matrix%b(1:PWGSMP__matrix%n)
      else
        PWGSMP__matrix%b => null()
      endif
    endif

    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__back_substitution() ..."

    return
  end subroutine PWGSMP__back_substitution
  
  
  
!  !> Apply iterative refinement
!  subroutine PWGSMP__iterative_refinement(rhs)
!    real*8, intent(inout) :: rhs(:)
!    integer :: i
!
!    include 'mpif.h'
!
!    if (PWSMP__verbose) WRITE(*,*) "Entering PWGSMP__iterative_refinement() ..."
!
!    if (PWGSMP__matrix%n > 0) then
!      if (PWSMP__copy_matrix) then
!        PWGSMP__matrix%b(1:PWGSMP__matrix%n) = rhs(1:PWGSMP__matrix%n)
!      else
!        ! associate pointer
!      endif
!    endif
!
!    PWGSMP__matrix%iparm(2) = 4
!    PWGSMP__matrix%iparm(3) = 4
!    call MPI_Barrier(MPI_COMM_WORLD, ierr)
!
!#ifdef USE_WSMP
!    call pwgsmp(PWGSMP__matrix%n,     PWGSMP__matrix%ir,    PWGSMP__matrix%jc,   PWGSMP__matrix%m,&
!                PWGSMP__matrix%b,     PWGSMP__matrix%ldb,   PWGSMP__matrix%nrhs, PWGSMP__matrix%rmsic,&
!                PWGSMP__matrix%iparm, PWGSMP__matrix%dparm)
!#endif
!    if (PWGSMP__matrix%iparm(64) .ne. 0) call wsmp_error_handler("PWGSMP back substitution failed.")
!
!    if (PWGSMP__matrix%n > 0) then
!      if (PWSMP__copy_matrix) then
!        rhs(1:PWGSMP__matrix%n) = PWGSMP__matrix%b(1:PWGSMP__matrix%n)
!      else
!        PWGSMP__matrix%b => null()
!      endif
!    endif
!
!    if (PWSMP__verbose) WRITE(*,*) "Exiting PWGSMP__iterative_refinement() ..."
!
!    return
!  end subroutine PWGSMP__iterative_refinement
  
  
  
  !> Print an error message and stop the code.
  subroutine wsmp_error_handler(errmsg)
    use mpi_mod
    character(len=*), intent(in) :: errmsg
    write(*,*) errmsg
    call MPI_Finalize(ierr)
    stop
  end subroutine wsmp_error_handler
  
end module wsmp_module

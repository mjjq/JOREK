module solve_mat_n

contains
#ifndef USE_PASTIX6
  ! -- For PaStiX solver before version 6.x
  subroutine pastix_bind_threads(my_id)
#else
  ! -- For PaStiX solver version 6.x
  subroutine pastix_bind_threads(my_id, thread_map)
#endif
    use pastix_module
    
    implicit none

#ifndef USE_PASTIX6
#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif
#endif

    integer, intent(in) :: my_id
#ifndef USE_PASTIX6
    ! -- For PaStiX solver before version 6.x
    integer*4, dimension(1:pastix_nthrd) :: thread_map
#else
    ! -- For PaStiX solver version 6.x
    integer(kind=c_int)    , dimension(1:pastix_nthrd), intent(out)   :: thread_map
#endif
    integer*4 k, packsize, procpernode

#if (defined(WORLDWAR2) || defined(USE_PASTIX6)) && defined(CORES_PER_NODE)
    procpernode = CORES_PER_NODE/pastix_nthrd
    packsize = CORES_PER_NODE/procpernode 
!    if (my_id .eq. 0) print *, "packsize", packsize, "procpernode", procpernode
    Do k = 1, pastix_nthrd
      thread_map(k) = mod(my_id * packsize,CORES_PER_NODE) + k-1
    end do
#ifndef USE_PASTIX6
    ! -- For PaStiX solver before version 6.x
    call pastix_fortran_bindthreads(pastix_data, pastix_nthrd, thread_map(1:))
#endif
#endif
  end subroutine pastix_bind_threads

  !> Solves the system of equation for each harmonic using mumps, pastix, or wsmp
  subroutine solve_matrix_n(my_id,i_tor,MPI_COMM_N,MPI_COMM_MASTER,solve_only)

    use tr_module
    use mod_parameters
    use mumps_module
    use wsmp_module
    use pastix_module
    use global_distributed_matrix
    use mpi_mod 
    use mod_clock
    use phys_module, only : index_now, use_BLR_compression, epsilon_BLR, just_in_time_BLR, pastix_blr_abs_tol
    use mod_coicsr
 
#ifdef USE_PASTIX6
    ! -- For PaStiX solver version 6.x
    use iso_c_binding
    use pastixf
    use pastix_enums
    use spmf
#endif
   
    implicit none

#ifndef USE_PASTIX6
    ! -- For PaStiX solver before version 6.x
#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif
#endif

#include "r3_info.h"

    integer, intent(in) :: my_id
    integer, dimension(:), intent(in) :: i_tor(:)
    integer, intent(in) :: MPI_COMM_N, MPI_COMM_MASTER
    logical, intent(in) :: solve_only

    integer :: i, j, k, my_id_n, n_cpu_n, ierr, my_id_master, n_cpu_master
    integer :: i_reduced, j_reduced, n_i, n_j, index, index1, index2
    type(clcktype) :: t_itstart, t0, t1, t2, t3
    real*8  :: tsecond
    real*8, allocatable :: RHS_tmp(:)
    !Split broadcast
    character*8 :: type
    real*8  :: DUMMY_REAL(1:1)
    integer :: DUMMY_INT (1:1)
    CHARACTER(LEN=128) :: fname
#ifdef USE_PASTIX6
    ! -- For PaStiX solver version 6.x
    integer(c_int)     :: pastix_info
    type(c_ptr)        :: pastix_rhs_ptr
    integer(kind=spm_int_t), dimension(:), pointer     :: pastix_colptr
    integer(kind=spm_int_t), dimension(:), pointer     :: pastix_rowptr
    real(kind=c_double)    , dimension(:), pointer     :: pastix_values
    integer(kind=c_int)    , dimension(1:pastix_nthrd) :: thread_map
#endif

    call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
    call tr_print_memsize("BeforeSolveN")
    call tr_debug_writei("smn_A_mumps_par%n",mumps_par%n)

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*      solve local matrix  (n)  *'
      write(*,*) my_id,'*********************************'

      if (use_mumps)  write(*,*) my_id,'*       using solver MUMPS      *'
      if (use_pastix) write(*,*) my_id,'*       using solver PastiX     *'
      if (use_wsmp)   write(*,*) my_id,'*       using solver WSMP       *'

      write(*,*) my_id,'*********************************'
    endif

    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)     ! the id of each cpu
    call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)     ! the number of cpus

    if (my_id_n .eq. 0) then
      call MPI_COMM_RANK(MPI_COMM_MASTER, my_id_master, ierr)     ! the id of each cpu
      call MPI_COMM_SIZE(MPI_COMM_MASTER, n_cpu_master, ierr)     ! the number of cpus
    endif


    if (.not. solve_only) then

      !---------------------------------------- column scaling 
      if (my_id_n .eq. 0) then

        if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
        call tr_allocate(column_scaling,1,mumps_par%N,"column_scaling",CAT_DMATRIX)

        column_scaling = 1.d-20
        do k=1,mumps_par%nz
          j = mumps_par%jcn(k)
          column_scaling(j) = min(max(column_scaling(j),abs(mumps_par%A(k))),1d20)
        enddo
        if (my_id .eq. 0) then
           write(fname,'(A,I6.6)')"column_scaling",index_now
           call tr_vdump(fname,column_scaling,mumps_par%N)
        end if
        !CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        write(*,'(2i4,A,2e12.4)') my_id,my_id_n,' COLUMN SCALING : ',minval(column_scaling),maxval(column_scaling)
        do k=1,mumps_par%nz
          j = mumps_par%jcn(k)
          mumps_par%A(k) = mumps_par%A(k) / column_scaling(j)
        enddo
      endif

      if (my_id_n .eq. 0) then                          ! elapsed time analysis start
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t0)
      endif


      if (use_mumps) then
#ifdef USE_MUMPS

        mumps_par%JOB = 1                                  ! Analysis, only needed when grid has changed

        mumps_par%icntl(7)  = mumps_ordering               ! ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS), default: 7
        mumps_par%icntl(8)  = 7                            ! row and column scaling
        mumps_par%icntl(14) = 30                           ! MAXS
        mumps_par%icntl(18) = 0

        if (use_BLR_compression) then
          mumps_par%icntl(35) = 1                          ! Block-low-rank (BLR) compression. 0: off (default), 1: automatic, 2: factorisation and solution, 3: only factorisation
          mumps_par%cntl(7)   = epsilon_BLR                ! Accuracy of BLR approximation
        endif

        call DMUMPS(mumps_par)

#endif
      else ! .not. use_mumps --> use_pastix or use_wsmp

        if (my_id_n .eq. 0) then           

          if (my_id_n .eq. 0) then                
            call MPI_Barrier(MPI_COMM_MASTER,ierr)
            call clck_time(t2)
          endif

#ifdef USE_BLOCK
          !---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
          !                             temporary solution before using blocks everywhere

          block_size  = n_var
          if (my_id .ne. 0) block_size = 2*n_var

          block_size2 = block_size**2
          !---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
          n_block   = mumps_par%n  / block_size
          nnz_block = mumps_par%nz / block_size2

          do i=1,nnz_block  
            mumps_par%irn(i) = (mumps_par%irn((i-1)*block_size2+1) - 1) / block_size + 1 
            mumps_par%jcn(i) = (mumps_par%jcn((i-1)*block_size2+1) - 1) / block_size + 1 
          enddo

          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(n_block+1))
          call coicsr2(n_block,nnz_block,mumps_par%A,mumps_par%IRN(1:nnz_block),mumps_par%JCN(1:nnz_block),block_size,sparskit_work)

          ! WARNING:  USE_BLOCK does not (yet) work with WSMP!!!
          if (use_wsmp) then
#ifdef USE_WSMP
            call PWGSMP__allocate(n_block, nnz_block, my_id_n)
            call PWGSMP__initialize_matrix(n_block, nnz_block,                                     &
              mumps_par%a, mumps_par%jcn, mumps_par%irn, my_id_n )
#endif
          endif

#else
          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(mumps_par%N + 1))

          call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

          if (use_wsmp) then
#ifdef USE_WSMP
            call PWGSMP__allocate(mumps_par%N, mumps_par%NZ, my_id_n)
            call PWGSMP__initialize_matrix(mumps_par%N, mumps_par%NZ,                              &
              mumps_par%a, mumps_par%jcn, mumps_par%irn, my_id_n )
#endif
          endif
#endif

          if (allocated(sparskit_work)) deallocate(sparskit_work)

          if (my_id_n .eq. 0) then
            call MPI_Barrier(MPI_COMM_MASTER,ierr) 
            call clck_time(t3)
            call clck_ldiff(t2,t3,tsecond)
            write(*,FMT_TIMING) my_id, '### Elapsed time coicsr :', tsecond
          endif

        else  ! (my_id_n > 0) below
#ifdef USE_WSMP
          if (use_wsmp) call PWGSMP__allocate(0, 0, my_id_n)
#endif
        endif ! end (my_id_n .eq. 0)



        ! --- Dstribute data to the MPI "slave" tasks (>0)
        !     (When using WSMP, this is *not necessary* in 0-master mode!)
        if ((.not. use_wsmp).and.(.not. pastix_smp_only)) then
!          call pastix_init_num_threads(my_id)

          call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
          call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
#ifdef USE_BLOCK
          call MPI_BCAST(block_size,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
          call MPI_BCAST(n_block,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
          call MPI_BCAST(nnz_block,1,MPI_INTEGER,0,MPI_COMM_N,ierr)
#endif
          if (my_id_n .gt. 0) then
            if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
            if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
            if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
            if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
            call tr_allocatep(mumps_par%irn,1,mumps_par%nz,"mumps_par%irn",CAT_DMATRIX)
            call tr_allocatep(mumps_par%jcn,1,mumps_par%nz,"mumps_par%jcn",CAT_DMATRIX)
            call tr_allocatep(mumps_par%a,1,mumps_par%nz,"mumps_par%a",CAT_DMATRIX)
            call tr_allocatep(mumps_par%rhs,1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)
          endif

          ! Split MPI_BCAST if MPI buffer beyond 2Go
          type='intIRN'
          call split_broadcast(type,MPI_COMM_N)
          type='intJCN'
          call split_broadcast(type,MPI_COMM_N)
          type='double'
          call split_broadcast(type,MPI_COMM_N)

#ifdef USE_PASTIX6
          ! -- For PaStiX solver version 6.x
          allocate(pastix_spm) ! Replace by tr_allocate etc.?!
          call spmInit(pastix_spm)

#ifdef USE_BLOCK
          pastix_spm%n           =  n_block
          pastix_spm%nnz         =  nnz_block
          pastix_spm%dof         =  block_size
#else
          pastix_spm%n           =  mumps_par%n
          pastix_spm%nnz         =  mumps_par%nz
          pastix_spm%dof         =  1
#endif
          call spmUpdateComputedFields(pastix_spm)
          call spmAlloc(pastix_spm)

          call c_f_pointer(pastix_spm%colptr,pastix_colptr, [pastix_spm%n+1])
          call c_f_pointer(pastix_spm%rowptr,pastix_rowptr, [pastix_spm%nnz])
          call c_f_pointer(pastix_spm%values,pastix_values, [mumps_par%nz])
              
          pastix_colptr      = mumps_par%jcn(1:pastix_spm%n+1)
          pastix_rowptr      = mumps_par%irn(1:pastix_spm%nnz)
          pastix_values      = mumps_par%A(1:mumps_par%nz)
#endif


        endif

        if  (.not. pastix_initialised)  then

          if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

            if (use_pastix) then

              call pastix_init_num_threads(my_id)

#ifndef USE_PASTIX6
              ! -- For PaStiX solver before version 6.x
              pastix_iparm(IPARM_MODIFY_PARAMETER) = API_NO         ! insert default values
              pastix_iparm(IPARM_START_TASK)       = API_TASK_INIT  ! initializse
              pastix_iparm(IPARM_END_TASK)         = API_TASK_INIT
#else
              ! -- For PaStiX solver version 6.x
              call pastixInitParam(pastix_iparm, pastix_dparm)
#endif

              if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%n,1,MPI_INTEGER,0,MPI_COMM_N,ierr)

#ifndef USE_PASTIX6
              ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
              call tr_allocate(pastix_perm_vars,1,n_block,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,1,n_block,"pastix_iperm_vars",CAT_UNKNOWN)

              call pastix_fortran(pastix_data,MPI_COMM_N,n_block,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
              call tr_allocate(pastix_perm_vars,1 ,mumps_par%n,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,1,mumps_par%n,"pastix_iperm_vars",CAT_UNKNOWN)

              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif
#endif

              ! pastix input parameters working in Pastix5 and Pastix6
              pastix_iparm(IPARM_VERBOSE)               = pastix_verb              
              pastix_iparm(IPARM_ITERMAX)               = pastix_iter                ! refinement : max number of iterations

              pastix_iparm(IPARM_FACTORIZATION)         = pastix_facto
              pastix_iparm(IPARM_THREAD_NBR)            = pastix_nthrd               ! number of threads
              pastix_iparm(IPARM_INCOMPLETE)            = pastix_ricar
              pastix_iparm(IPARM_LEVEL_OF_FILL)         = pastix_iluk
              pastix_dparm(DPARM_EPSILON_REFINEMENT)    = pastix_epsilon             ! error level refinement
              pastix_dparm(DPARM_EPSILON_MAGN_CTRL)     = pastix_pivot               ! pivot threshold
#ifdef USE_BLOCK
              pastix_iparm(IPARM_DOF_NBR)               = block_size                 ! block size
#else
              pastix_iparm(IPARM_DOF_NBR)               = 1
#endif



#ifndef USE_PASTIX6
              ! -- For PaStiX solver before version 6.x
              pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
              pastix_iparm(IPARM_SYM)                   = pastix_sym
              pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg
#ifdef WORLDWAR2
#ifdef FUNNELED
              pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
              pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif
#endif


#else
              ! -- For PaStiX solver version 6.x
              pastix_iparm(IPARM_MTX_TYPE)              = pastix_sym
              pastix_iparm(IPARM_AMALGAMATION_LVLCBLK)  = pastix_amalg

! TEMPORARY: not yet relevant for Pastix6 as MPI parallelisation is not implemented
!#ifdef FUNNELED
!              pastix_iparm(IPARM_THREAD_COMM_MODE)      = PastixThreadFunneled
!#else
!              pastix_iparm(IPARM_THREAD_COMM_MODE)      = PastixThreadMultiple
!#endif
              ! BLR Compression
              if (use_BLR_compression) then
                if (just_in_time_BLR) then
                  pastix_iparm(IPARM_COMPRESS_WHEN)     = PastixCompressWhenEnd ! Just-in-Time (speed optimal)
                else 
                  pastix_iparm(IPARM_COMPRESS_WHEN)     = PastixCompressWhenBegin ! Minimal-memory (default)
                endif
                if (pastix_blr_abs_tol) then
                  pastix_iparm(IPARM_COMPRESS_RELTOL)     = 0
                else
                  pastix_iparm(IPARM_COMPRESS_RELTOL)     = 1
                end if
                pastix_dparm(DPARM_COMPRESS_TOLERANCE)  = epsilon_BLR

!!               Additional PaStiX compression parameters (currently set to their default values)
!                pastix_iparm(IPARM_COMPRESS_ORTHO)      = PastixCompressOrthoCGS
!                pastix_iparm(IPARM_COMPRESS_METHOD)     = PastixCompressMethodPQRCP
!                pastix_iparm(IPARM_COMPRESS_MIN_WIDTH)  = 120
!                pastix_iparm(IPARM_COMPRESS_MIN_HEIGHT) = 20
!                pastix_dparm(DPARM_COMPRESS_MIN_RATIO)  = 1.0
              endif

              ! initialise PaStiX (and bind threads if desired, else automatic binding)
#ifdef CORES_PER_NODE
              call pastix_bind_threads(my_id, thread_map)
              call pastixInitWithAffinity(pastix_data, 0, pastix_iparm, pastix_dparm, thread_map)    ! TEMPORARY: 0 should be pastix_comm but pastix6 is not yet MPI parallelised!
#else
              call pastixInit(pastix_data, 0, pastix_iparm, pastix_dparm)    ! TEMPORARY: 0 should be pastix_comm but pastix6 is not yet MPI parallelised!
#endif

#endif

             else if (use_wsmp) then
#ifdef USE_WSMP
              call PWGSMP__initialize_solver(my_id_n, MPI_COMM_N)
#endif
            end if

           pastix_initialised = .true.

          endif

        endif !.not. pastix_initialised



        if (.not. pastix_analysed) then
          if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

            if (use_pastix) then
#ifndef USE_PASTIX6
              ! -- For PaStiX solver before version 6.x
              pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
              pastix_iparm(IPARM_START_TASK) = API_TASK_ORDERING
              pastix_iparm(IPARM_END_TASK)   = API_TASK_ANALYSE
!              pastix_iparm(IPARM_BINDTHRD)   = API_NO
#ifdef USE_BLOCK
              call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
                mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif

#else
              ! -- For PaStiX solver version 6.x
#ifdef USE_BLOCK
! ############################################################################
! ####### these lines can be replaced by pastix_task_analyze in the future,
! ####### as soon as PaStiX 6 supports multiple dofs in all solver steps.
! ############################################################################
              call pastix_subtask_order(pastix_data,pastix_spm,pastix_myorder,pastix_info)
              call pastix_subtask_symbfact(pastix_data,pastix_info)
              call pastix_subtask_reordering(pastix_data,pastix_info)

              ! Expand spm matrix and pastix analysis substructures because rest of Pastix6 cannot handle multiple dofs (yet)
              call pastixExpand(pastix_data,pastix_spm)
             
              call pastix_subtask_blend(pastix_data,pastix_info)
! ############################################################################
! ####### end these lines can be replaced...
! ############################################################################
#else
              call pastix_task_analyze(pastix_data,pastix_spm,pastix_info)
#endif
#endif
            else if (use_wsmp) then
              ! do nothing
            endif

            pastix_analysed = .true.
#if (defined(USE_PASTIX6) && defined(USE_BLOCK))
            pastix_analysed = .false. ! Necessary for now such that the spm expansion is done every time step. 
                                      ! Can be removed once the PaStiX team has implemented multi-dof for all pastix_subtasks.
#endif
          endif
        endif ! .not. pastix_analysed
      endif   ! (else, use_mumps)


      if (my_id_n .eq.0) then                            ! elapsed time analysis end
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t1)
         call clck_ldiff(t0,t1,tsecond)
         write(*, FMT_TIMING) my_id,' ## Elapsed time, analysis :',tsecond
         call clck_time(t0)                              ! elapsed time facto start 
      endif

      if (use_mumps) then
#ifdef USE_MUMPS

        mumps_par%JOB = 2                                   ! factorisation

        call DMUMPS(mumps_par)

        if (my_id_n .eq.0)   write(*,'(i3,A,i8)')    my_id,' MUMPS, mem       : ',mumps_par%info(16)
#endif
      elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

        if (use_pastix) then

#ifndef USE_PASTIX6
          ! -- For PaStiX solver before version 6.x
          pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
          pastix_iparm(IPARM_START_TASK) = API_TASK_NUMFACT
          pastix_iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
#if defined(WORLDWAR2) && defined(CORES_PER_NODE)
          pastix_iparm(IPARM_BINDTHRD)   = API_BIND_TAB
#endif
          call pastix_bind_threads(my_id)
#ifdef USE_BLOCK
          call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
            mumps_par%jcn, mumps_par%irn, mumps_par%A, &
            pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#else	   
          call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
            pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif

#else
          ! -- For PaStiX solver version 6.x
          call pastix_task_numfact(pastix_data,pastix_spm,pastix_info)

          call spmExit(pastix_spm)
          deallocate(pastix_spm)
#endif

        else if (use_wsmp) then
#ifdef USE_WSMP
          call PWGSMP__LU_factorization(my_id_n)
#endif
        end if

      endif
     
      if (my_id_n .eq.0) then                            ! elapsed time facto end
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t1)
         call clck_ldiff(t0,t1,tsecond)
         write(*, FMT_TIMING) my_id,' ## Elapsed time, facto :',tsecond
         call clck_time(t0)
      end if
   endif
   call tr_debug_writei("smn_B_mumps_par%n",mumps_par%n)


   if (my_id_n .eq. 0) then                          ! elapsed time solve start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
      call clck_time(t0)
   endif
   if (use_mumps) then
#ifdef USE_MUMPS
      mumps_par%JOB = 3                                   ! Solve

      call DMUMPS(mumps_par)

#endif
    elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

      if (use_pastix) then
        if (.not. pastix_smp_only) then
           call tr_debug_writei("smn_C_mumps_par%n",mumps_par%n)
           call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
        end if
        if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
        if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
        if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)

        call tr_locvnorms("smn_rhs",mumps_par%rhs,mumps_par%n)

#ifndef USE_PASTIX6
        ! -- For PaStiX solver before version 6.x
        pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
        pastix_iparm(IPARM_START_TASK) = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK)   = pastix_endsolve
!        pastix_iparm(IPARM_BINDTHRD)   = API_NO
#ifdef USE_BLOCK
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                &
!             mumps_par%jcn,mumps_par%irn,mumps_par%A, &
             DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
             pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,&
!             mumps_par%jcn,mumps_par%irn,mumps_par%A, &
          DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#endif

#else
       ! -- For PaStiX solver version 6.x
       pastix_rhs_ptr = c_loc(mumps_par%rhs)
#ifdef USE_BLOCK
       call pastix_task_solve(pastix_data,1,pastix_rhs_ptr,n_block,pastix_info)
#else
       call pastix_task_solve(pastix_data,1,pastix_rhs_ptr,mumps_par%n,pastix_info)
#endif
#endif

      else if (use_wsmp) then
#ifdef USE_WSMP
        call PWGSMP__back_substitution(mumps_par%rhs, my_id_n)
#endif
     end if
   endif

    if (my_id_n .eq.0) then                            ! elapsed time solve end
       call MPI_Barrier(MPI_COMM_MASTER,ierr)
       call clck_time(t1)
       call clck_ldiff(t0,t1,tsecond)
       write(*, FMT_TIMING) my_id,' ## Elapsed time, solve :',tsecond
       call clck_time(t0)
    end if


    if (my_id_n .eq. 0) then

      !------------------------------------------ undo column scaling
      do k=1,mumps_par%n
        mumps_par%rhs(k) =  mumps_par%rhs(k) / column_scaling(k)
      enddo

      if (allocated(deltas)) call tr_deallocate(deltas,"deltas",CAT_PRECOND)
      call tr_allocate(deltas,1,ndof_glob,"deltas",CAT_PRECOND)
      deltas = 0.d0

      call tr_allocate(rhs_tmp,1,ndof_glob,"rhs_tmp",CAT_PRECOND)

      rhs_tmp = 0.d0

      if (my_id .eq. 0 ) then
        !        rhs_tmp(1:ndof_glob:n_tor) = mumps_par%rhs(1:mumps_par%n)
        do i=0, mumps_par%n-1
          rhs_tmp(1+i*n_tor)=mumps_par%rhs(1+i)
        end do
      else
        !        rhs_tmp(2*i_tor(my_id+1)-2:ndof_glob:n_tor) = mumps_par%rhs(1:mumps_par%n:2)
        !        rhs_tmp(2*i_tor(my_id+1)-1:ndof_glob:n_tor) = mumps_par%rhs(2:mumps_par%n:2)
        do i=0, mumps_par%n/2-1
          rhs_tmp(2*i_tor(my_id+1)-2+i*n_tor) = mumps_par%rhs(1+i*2)
          rhs_tmp(2*i_tor(my_id+1)-1+i*n_tor) = mumps_par%rhs(2+i*2)
        end do

      endif

      call MPI_AllReduce(RHS_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
      call tr_deallocate(rhs_tmp,"rhs_tmp",CAT_PRECOND)

      call tr_locvnorms("smn_res",mumps_par%rhs,mumps_par%n)
      call tr_locvnorms("smn_delta",deltas,ndof_glob)
    endif
#ifndef USE_PASTIX6
    ! -- For PaStiX solver before version 6.x
    call tr_set_precondmem(pastix_dparm(DPARM_MEM_MAX)) ! DPARM_MEM_MAX DEPRECATED IN PASTIX6: how to change this?
    ! ############### This should be looked at at some point
#endif
    call tr_print_memsize("AfterSolveN")
    call r3_info_end (r3_info_index_0)         ! timing
    return
  end subroutine solve_matrix_n
end module solve_mat_n

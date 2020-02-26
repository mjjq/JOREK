!> Solve step of the local matrices for each toroidal harmonic (preconditioner for gmres)
subroutine gmres_precondition(x,y,i_tor,my_id,my_id_n,MPI_COMM_MASTER,MPI_COMM_N)

  use tr_module 
  use mod_parameters
  use mumps_module
  use pastix_module
  use wsmp_module
  use global_distributed_matrix
  use mpi_mod
  use mod_clock
  use phys_module, only: use_pastix, use_mumps, use_strumpack

#ifdef USE_PASTIX6
! -- For PaStiX solver version 6.x
  use pastixf
  use pastix_enums
  use spmf
#endif

#ifndef USE_PASTIX6
! -- For PaStiX solver before version 6.x
#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif
#endif

#ifdef USE_STRUMPACK
  use spk_module
#endif

  implicit none

  integer             :: my_id, my_id_n, MPI_COMM_MASTER, MPI_COMM_N, ierr, i, k, i_tor(*), n_dof
  integer             :: my_id_master
  real*8              :: x(*), y(*)
  real*8, allocatable :: y_tmp(:), Rsnd_buffer(:)
  integer             :: index_snd, n_loc_n, n_cpu, n_cpu_n, M_cpu, ifactor, in, j, idisp, index_rcv
  integer, allocatable :: send_counts(:), send_disp(:), recv_counts(:), recv_disp(:)
  type(clcktype)       :: t0, t1
  real*8               :: tsecond
  real*8, allocatable :: buffer(:)
  integer             :: ibuf_size, status(MPI_STATUS_SIZE)
  
  real*8  :: DUMMY_REAL(1:1)
  integer :: DUMMY_INT (1:1)

#ifdef USE_PASTIX6
! -- For PaStiX solver version 6.x
  integer(c_int)     :: pastix_info
  type(c_ptr)        :: pastix_rhs_ptr
#endif

!write(*,*) my_id,my_id_n,' GMRES preconditioning ',MPI_COMM_WORLD

  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)     ! the number of cpus

  if (my_id_n .eq. 0) call MPI_COMM_RANK(MPI_COMM_MASTER, my_id_master, ierr)      ! the id of each cpu in the MASTER group

  call clck_time(t0)

  n_dof    = ndof_glob

  n_loc_n  = n_dof / n_tor
  M_cpu    = n_cpu / ((n_tor+1)/2)

  ifactor = 2
  if (my_id < M_cpu) ifactor = 1

  if (my_id .eq. 0) then

    call tr_allocate(Rsnd_buffer,1,n_dof,"Rsnd_buffer",CAT_GMRES,.false.)
    call tr_allocate(send_counts,1,n_cpu/M_cpu,"send_counts",CAT_GMRES)
    call tr_allocate(send_disp,1,n_cpu/M_cpu,"send_disp",CAT_GMRES)
    Rsnd_buffer(1:n_loc_n) = x(1:n_dof:n_tor)

    do in=2, (n_tor+1)/2

      index_snd = n_loc_n + (in-2)*2*n_loc_n

      Rsnd_buffer(index_snd+1:index_snd+2*n_loc_n:2) = x(2*in-2:n_dof:n_tor)
      Rsnd_buffer(index_snd+2:index_snd+2*n_loc_n:2) = x(2*in-1:n_dof:n_tor)

    enddo

    send_counts(1) = n_loc_n
    send_disp(1)   = 0
    idisp          = send_counts(1)

    do j=2,n_cpu/M_cpu

      send_counts(j) = 2*n_loc_n
      send_disp(j)   = idisp
      idisp          = idisp + send_counts(j)

    enddo
  else
    call tr_allocate(Rsnd_buffer,1,1,"Rsnd_buffer",CAT_GMRES)
    call tr_allocate(send_counts,1,1,"send_counts",CAT_GMRES)
    call tr_allocate(send_disp,1,1,"send_disp",CAT_GMRES)
  endif

  if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)

  if (my_id_n .eq. 0) then

    call tr_allocatep(mumps_par%rhs,1,ifactor*n_loc_n,"mumps_par%rhs",CAT_DMATRIX)

    call mpi_scatterv(Rsnd_buffer,send_counts,send_disp,MPI_DOUBLE_PRECISION, &
                    mumps_par%rhs,ifactor*n_loc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_MASTER,ierr)

!------------- mpi_scatterv alternative
!  ibuf_size = 8*n_dof  
!  allocate(buffer(ibuf_size))
!  call mpi_buffer_attach(buffer,ibuf_size,ierr)
!  if (my_id_master .eq. 0) mumps_par%rhs(1:n_loc_n) = Rsnd_buffer(1:n_loc_n)
!  do i=2,(n_tor+1)/2 
!    if (my_id_master .eq. 0) then    
!      idisp = n_loc_n + 1  + (i-2)*2*n_loc_n 
!      call mpi_bsend(Rsnd_buffer(idisp),2*n_loc_n,MPI_DOUBLE_PRECISION,i-1,i-1,MPI_COMM_MASTER,ierr)    
!    endif      
!    if (my_id_master .eq. i-1) then
!      call mpi_recv(mumps_par%rhs,2*n_loc_n,MPI_DOUBLE_PRECISION,0,i-1,MPI_COMM_MASTER,status,ierr)      
!    endif    
!  enddo
!  call mpi_buffer_detach(buffer,ibuf_size,ierr)
!------------------- end alternative

  endif

  call tr_deallocate(Rsnd_buffer,"Rsnd_buffer",CAT_GMRES)
  call tr_deallocate(send_counts,"send_counts",CAT_GMRES)
  call tr_deallocate(send_disp,"send_disp",CAT_GMRES)


!call clck_time(t1)
!call clck_ldiff(t0,t1,tsecond)
!t0 = t1
!if (my_id_n .eq. 0)  then
!   write(*,FMT_TIMING) my_id, '## Elapsed time precondition1 :', tsecond
!end if

#ifdef USE_MUMPS
  if (use_mumps) then
    mumps_par%JOB = 3                                   ! Solve
    call DMUMPS(mumps_par)
  endif
#endif

#if defined(USE_PASTIX)||defined(USE_PASTIX6)
  if (use_pastix) then
    if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then
  
      if (.not. associated(mumps_par%rhs)) then
        !    write(*,*) ' gmres: RHS not allocated!',my_id, my_id_n
        call tr_allocatep(mumps_par%rhs,1,ifactor*n_loc_n,"mumps_par%rhs",CAT_DMATRIX)
      endif
     
      if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%rhs,ifactor*n_loc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
  
        ! pastix input parameters working in Pastix5 and Pastix6
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
        pastix_iparm(IPARM_START_TASK)            = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK)              = pastix_endsolve
        pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
        pastix_iparm(IPARM_SYM)                   = pastix_sym
        pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg
  
#else
        ! -- For PaStiX solver version 6.x
        pastix_iparm(IPARM_MTX_TYPE)              = pastix_sym
        pastix_iparm(IPARM_AMALGAMATION_LVLCBLK)  = pastix_amalg
#endif
       
#ifndef USE_PASTIX6
        ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                        &
             !mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                   DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else      
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n, DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
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
    endif
  endif ! use_pastix
#endif /* defined(USE_PASTIX)||defined(USE_PASTIX6) */

#ifdef USE_WSMP
  if (use_wsmp) then
    call PWGSMP__back_substitution(mumps_par%rhs, my_id_n)
  endif  
#endif

#ifdef USE_STRUMPACK
  if (use_strumpack) then
    if (.not. associated(mumps_par%rhs)) then
       call tr_allocatep(mumps_par%rhs,1,ifactor*n_loc_n,"mumps_par%rhs",CAT_DMATRIX)
    endif
   
    call MPI_BCAST(mumps_par%rhs,ifactor*n_loc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)      
    call f2spk(ifactor*n_loc_n,mumps_par%nz,mumps_par%irn,mumps_par%jcn,mumps_par%A,mumps_par%rhs,MPI_COMM_N,3)
  endif  
#endif 

!call clck_time(t1)
!call clck_ldiff(t0,t1,tsecond)
!t0 = t1
!if (my_id_n .eq. 0)  then
!   write(*,FMT_TIMING) my_id, '## Elapsed time precondition2 :', tsecond
!end if

!write(*,'(2i3,A,2e16.8)') my_id,my_id_n,' precond : rhs after : ',maxval(mumps_par%rhs),minval(mumps_par%rhs)


  if (my_id_n .eq. 0) then
  
    if (.not.use_strumpack) then
  !------------------------------------------ undo column scaling
      do k=1,mumps_par%n
        mumps_par%rhs(k) =  mumps_par%rhs(k) / column_scaling(k)
      enddo
    endif
  
    call tr_allocate(y_tmp,1,n_dof,"y_tmp",CAT_GMRES,.false.)
    call tr_allocate(recv_counts,1,n_cpu/M_cpu,"recv_counts",CAT_GMRES)
    call tr_allocate(recv_disp,1,n_cpu/M_cpu,"recv_disp",CAT_GMRES)
  
    y_tmp(1:n_dof) = 0.d0
  
    recv_counts(1) = n_loc_n
    do i=2,(n_tor+1)/2
      recv_counts(i) = 2*n_loc_n
    enddo
  
    recv_disp(1) = 0
    do i=2,(n_tor+1)/2
      recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
    enddo
  
    call mpi_gatherv(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION, &
                     y_tmp,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,0,MPI_COMM_MASTER,ierr)
  
  !----------------------------- mpi_gatherv alternative
  !  if (my_id_master .eq. 0)  y_tmp(1:n_loc_n) = mumps_par%rhs(1:n_loc_n) 
  !  call mpi_buffer_attach(buffer,ibuf_size,ierr)
  !  do i=2,(n_tor+1)/2 
  !    if (my_id_master .eq. i-1) then
  !      call mpi_bsend(mumps_par%rhs,2*n_loc_n,MPI_DOUBLE_PRECISION,0,i-1,MPI_COMM_MASTER,ierr)
  !    endif  
  !    if (my_id_master .eq. 0) then    
  !      idisp = n_loc_n + 1  + (i-2)*2*n_loc_n 
  !      call mpi_recv(y_tmp(idisp),2*n_loc_n,MPI_DOUBLE_PRECISION,i-1,i-1,MPI_COMM_MASTER,status,ierr)
  !    endif    
  !  enddo
  !  call mpi_buffer_detach(buffer,ibuf_size,ierr)
  !  deallocate(buffer)
  !--------------------------- end alternative
  
    if (my_id .eq. 0) then
  
      y(1:n_dof:n_tor) = y_tmp(1:n_loc_n)
  
      do in=2, (n_tor+1)/2
  
        index_rcv = n_loc_n + (in-2)*2*n_loc_n
  
        y(2*in-2:n_dof:n_tor) = y_tmp(index_rcv+1:index_rcv+2*n_loc_n:2)
        y(2*in-1:n_dof:n_tor) = y_tmp(index_rcv+2:index_rcv+2*n_loc_n:2)
  
      enddo
  
    endif
  
    call tr_deallocate(y_tmp,"y_tmp",CAT_GMRES)
    call tr_deallocate(recv_counts,"recv_counts",CAT_GMRES)
    call tr_deallocate(recv_disp,"recv_disp",CAT_GMRES)
  endif

!call clck_time(t1)
!call clck_ldiff(t0,t1,tsecond)
!t0 = t1
!if (my_id_n .eq. 0)  then
!   write(*,FMT_TIMING) my_id, '## Elapsed time precondition3 :', tsecond
!end if
  return
end subroutine gmres_precondition

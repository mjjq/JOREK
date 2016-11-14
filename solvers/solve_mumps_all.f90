subroutine solve_mumps_all(my_id)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using mumps with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use tr_module 
use mumps_module
use global_distributed_matrix
use mpi_mod
use mod_clock

implicit none

real*8,allocatable       :: column_local(:)
real*8                   :: tsecond, t_analysis_0, t_analysis_1, t_fact_0, t_fact_1
type(clcktype)           :: t0, t1
integer                  :: k, j, ierr, my_id
#ifdef USE_MUMPS
!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*      solve global matrix      *'
!write(*,*) my_id,'*********************************'

mumps_par%A_loc   => A_glob(1:nz_glob)
mumps_par%irn_loc => irn_glob(1:nz_glob)
mumps_par%jcn_loc => jcn_glob(1:nz_glob)
mumps_par%rhs     => rhs_glob(1:nz_glob)

mumps_par%n      = ndof_glob
mumps_par%nz_loc = nz_glob

if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
if (allocated(column_local))    call tr_deallocate(column_local,"column_local",CAT_DMATRIX)
call tr_allocate(column_scaling,1,mumps_par%N,"column_scaling",CAT_DMATRIX)
call tr_allocate(column_local,1,mumps_par%N,"column_local",CAT_DMATRIX)

column_local = 1.d-20;   column_scaling = 1.d-20
do k=1,mumps_par%nz_loc
  j = mumps_par%jcn_loc(k)
  column_local(j) = max(column_local(j),abs(mumps_par%A_loc(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
do k=1,mumps_par%nz_loc
  j = mumps_par%jcn_loc(k)
  mumps_par%A_loc(k) = mumps_par%A_loc(k) / column_scaling(j)
enddo

mumps_par%JOB = 1                                  ! Analysis, only needed when grid has changed
mumps_par%icntl(7)  = 7
mumps_par%icntl(8)  = 7                            ! row and column scaling  7: automatic scaling
mumps_par%icntl(18) = 3
mumps_par%icntl(14) = 50                           ! MAXS

call clck_time(t0)

call DMUMPS(mumps_par)

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  then
   write(*,FMT_TIMING) my_id, '## Elapsed time mumps analyis :', tsecond
end if

call clck_time(t0)
mumps_par%JOB = 2                                   ! factorisation
call DMUMPS(mumps_par)

mumps_par%JOB = 3                                   ! Solve
call clck_time(t1)
call DMUMPS(mumps_par)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  then
   write(*,FMT_TIMING) my_id, '## Elapsed time mumps fact/solve :', tsecond
end if

do k=1,mumps_par%n
  deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
!  write(*,*) k,deltas(k)
enddo
#endif
return
end

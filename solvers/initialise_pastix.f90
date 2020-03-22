subroutine initialise_pastix
!***********************************************************************
! initialise the pastix library                                        *
!***********************************************************************
use pastix_module
use mumps_module
use phys_module, only: use_pastix

implicit none

if (.not. use_pastix) return

write(*,*) '***********************************'
write(*,*) '* initialise PastiX               *'
write(*,*) '***********************************'

pastix_iparm(1)     = 0                                     ! insert default values
pastix_iparm(2)     = 0                                     ! initializse
pastix_iparm(3)     = 0

pastix_iparm(31) = pastix_facto
pastix_iparm(35) = pastix_nthrd                  ! thread/mpi
pastix_iparm(39) = pastix_rhs
pastix_iparm(41) = pastix_sym

pastix_iparm(42) = pastix_ricar
pastix_iparm(37) = pastix_iluk
pastix_iparm(14) = pastix_amalg

!call my_pastix(MPI_COMM_N,mumps_par%n,mumps_par%irn,mumps_par%jcn,mumps_par%A,   &
!                   pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

pastix_iparm(6)  = pastix_iter
pastix_dparm(6)  = pastix_epsilon
pastix_dparm(11) = pastix_pivot             ! pivot threshold?

return
end

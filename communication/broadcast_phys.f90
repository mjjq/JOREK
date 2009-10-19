subroutine Broadcast_phys(my_id)
!----------------------------------------------------------
! subroutine to broadcast all the nodes in the point_list
!----------------------------------------------------------
use phys_module

implicit none

include 'mpif.h'               ! MPI fortran include file
integer              :: my_id, ierr, INT_EXT, IDBL_EXT, ILOG_EXT,position, bufsize
integer, allocatable :: buffer(:)

!----------------------------------- one line would be enough if only MPI_TYPE_STRUCT would work on IXIA
!call MPI_BCAST(phys_list,1,MPI_phys,0,MPI_COMM_WORLD,ierr)


call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,IDBL_EXT,ierr)
call MPI_TYPE_EXTENT(MPI_INTEGER,INT_EXT,ierr)
call MPI_TYPE_EXTENT(MPI_LOGICAL,ILOG_EXT,ierr)

bufsize = ( 132 * IDBL_EXT + (4+n_tor) * INT_EXT + 5 * ILOG_EXT )

if (allocated(buffer)) deallocate(buffer)
allocate(buffer(bufsize/ INT_EXT))

if (my_id .eq. 0) then
  position = 0

  call MPI_PACK(eta,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(visco,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(visco_par,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(tstep,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(F0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(GAMMA,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(Q_bar,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(sigma,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  call MPI_PACK(zjz_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(zjz_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(zj_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(T_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(T_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(T_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(Ti_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(Ti_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(Ti_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(Te_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(Te_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(T_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(rho_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(rho_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(rho_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(FF_0,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(FF_1,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(FF_coef,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(heatsource,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(heatsource_i,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(heatsource_e,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(particlesource,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
  call MPI_PACK(ZK_perp,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(ZK_par,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(ZK_i_perp,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(K_i_par,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(ZK_e_perp,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(K_e_par,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(D_perp,10,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(D_par,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(eta_num,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(visco_num,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(visco_par_num,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(D_perp_num,1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(nstep,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(n_flux,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(mode,n_tor,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(index_start,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(index_now,1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

  call MPI_PACK(restart,1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(regrid,1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(import_equil,1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(xpoint,1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  call MPI_PACK(freeboundary,1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

  position = 0

  call MPI_UNPACK(buffer,bufsize,position,eta,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_par,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,tstep,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,F0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,GAMMA,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Q_bar,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,sigma,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,zjz_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zjz_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,zj_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,T_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,T_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,T_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,Ti_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Ti_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Ti_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,Te_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Te_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,Te_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,rho_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rho_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,rho_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,FF_0,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,FF_1,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,FF_coef,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,heatsource,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,heatsource_i,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,heatsource_e,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,particlesource,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,ZK_perp,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ZK_par,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ZK_i_perp,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,K_i_par,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,ZK_e_perp,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,K_e_par,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_perp,10,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_par,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,eta_num,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_num,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,visco_par_num,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,D_perp_num,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,nstep,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,n_flux,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,mode,n_tor,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,index_start,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,index_now,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

  call MPI_UNPACK(buffer,bufsize,position,restart,1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,regrid,1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,import_equil,1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,xpoint,1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  call MPI_UNPACK(buffer,bufsize,position,freeboundary,1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)

endif

if (allocated(buffer)) deallocate(buffer)

return
end

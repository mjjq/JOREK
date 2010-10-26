subroutine Poisson(my_id,itype,node_list,element_list,ivar_in,ivar_out,i_harm,xpoint)
!---------------------------------------------------------------
! collect the element matrices into one large sparse matrix
! in coordinate format
!---------------------------------------------------------------
use data_structure
use mumps_module
use pastix_module
implicit none
include 'mpif.h'

interface
   subroutine find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis)
     !-----------------------------------------------------------------------
     !-----------------------------------------------------------------------
     use data_structure

     type (type_node_list)    :: node_list
     type (type_element_list) :: element_list

     real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
     integer :: i_elm_axis
   end subroutine find_axis
end interface

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
type (type_element)      :: element_father
type (type_node)         :: nodes_father(n_vertex_max)
real*8   :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1))
real*8   :: zbig, Z_xpoint, psi_axis, psi_bnd, psi_xpoint, R_xpoint, s_xpoint, t_xpoint
real*8   :: R_axis, Z_axis, s_axis, t_axis, amix
integer  :: ierr, my_id, itype, ivar_in, ivar_out, i_harm, i_elm_axis, i_elm_xpoint
integer  :: n_AA, nz_AA, nz_AA_old, n_border, ilarge, ife, iv, i,j,k,l
integer  :: n_elements, inode, index_large_i, knode, index_large_k, index_ij, index_kl, index, index_i, check_data, nnz
logical  :: xpoint
real*8, dimension(4,4)	 :: H, H_s, H_t, H_st
real*8, dimension(2,4) 	 :: c, dc_ds, dc_dt, d2c_dsdt					   
real*8			 :: lambda, mu	
real*8			 :: Psi,dPsi_ds,dPsi_dt,d2Psi_dsdt
real*8			 :: dX_ds, dX_dt, dY_ds, dY_dt, d2X_dsdt, d2Y_dsdt, h_u, h_v, h_w 
integer			 :: i_father, inode_father,Index_elm    
integer, dimension(n_vertex_max)  :: pr
integer, dimension(2)		  :: parent
integer, dimension(n_vertex_max) ::  node_out
!integer, dimension(node_list%n_nodes) :: active_node
integer  ::n_active_nodes
if (my_id .eq. 0) then

  write(*,*) '**************************************'
  write(*,*) '*            Poisson                 *'
  write(*,*) '**************************************'
  write(*,*) ' i_type     : ',itype
  write(*,*) ' n_elements : ',element_list%n_elements
  write(*,*) ' n_nodes    : ',node_list%n_nodes

  ! n_active_nodes = 0
  
 ! do i = 1, node_list%n_nodes
       
      
           ! if((node_list%node(i)%constrained==.false.) )  then
                ! n_active_nodes = n_active_nodes + 1
	        ! active_node(i) = n_active_nodes	! Correspondance entre numéro du noeud et position 
		 					! dans la liste des noeuds actifs
           !do k=1,n_order+1
            !node_list%node(i)%index(k) =  (n_order+1)*(active_node(i)-1)+k 
           !print*,"idex",i,node_list%node(i)%index(k)
          ! enddo	
         !else
	    !  active_node(i) = 0							                          
       !end if 
      
  !end do
 !stop
  
  nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2

  n_border = 0
  do i=1,node_list%n_nodes
   if (node_list%node(i)%boundary .eq. 1) n_border = n_border+2
   if (node_list%node(i)%boundary .eq. 2) n_border = n_border+2
   if (node_list%node(i)%boundary .eq. 3) n_border = n_border+3
  enddo

  nz_AA = nz_AA + n_border

  n_AA  = node_list%n_nodes * (n_order+1)

  n_AA = 0
  do inode = 1, node_list%n_nodes
    n_AA = max(n_AA,node_list%node(inode)%index(4))
  enddo

  write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)
  write(*,*) ' number of boundary nodes: ',n_border

  if (.not. associated(mumps_par%A))     allocate(mumps_par%A(nz_AA))
  if (.not. associated(mumps_par%rhs))   allocate(mumps_par%rhs(n_AA))
  if (.not. associated(mumps_par%irn))   allocate(mumps_par%irn(nz_AA))
  if (.not. associated(mumps_par%jcn))   allocate(mumps_par%jcn(nz_AA))

  mumps_par%irn = 0
  mumps_par%jcn = 0
  mumps_par%A   = 0.d0
  mumps_par%RHS = 0.d0

  n_elements = element_list%n_elements

  ilarge=0

  if (itype .lt. 0) then
     i_elm_axis = 1 ! XL : uninitilised value... So let's say it'll be 1, to look in the first case of element array in interp.f90...
     call find_axis(node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis)

    psi_bnd = 0.d0
    if (xpoint) then
      call find_xpoint(node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint)
      psi_bnd = psi_xpoint
    endif
  endif

  amix = 0.d0
  if (itype .eq. -1) amix = 0.75d0

  do ife =1, n_elements

    element = element_list%element(ife)   
    i_father= element_list%element(ife)%father

    if( i_father.ne.0) then

      element_father = element_list%element(i_father)

    endif

    do iv = 1, n_vertex_max
      
      if( i_father.ne.0) then
         inode_father=element_father%vertex(iv)
         nodes_father(iv) = node_list%node(inode_father)
      endif

      inode     = element%vertex(iv)
      nodes(iv) = node_list%node(inode)

    enddo

    if (itype .eq. -1) then

      call element_matrix_GS(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

    elseif (itype .eq. -2) then

      call element_matrix_GS_inverse(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

    else

      call element_matrix_Poisson(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

    endif!n_active_nodes = 0
  
  !do i = 1, node_list%n_nodes
       
      
            !if((node_list%node(i)%constrained==.false.) )  then
                ! n_active_nodes = n_active_nodes + 1
	        ! active_node(i) = n_active_nodes	
           !do k=1,n_order+1
            !node_list%node(i)%index(k) =  (n_order+1)*(active_node(i)-1)+k 
           
           !enddo	
        ! else
	     ! active_node(i) = 0							                          
       !end if 
  
  !end do
!return
  !write(*,*) 'Noeuds actifs = ', n_active_nodes,node_list%n_nodes
  
   ! if(refinement ==.true.) then

       call Chgmt_node(ife,element,nodes,element_father,nodes_father,ELM,RHS,node_out) ! Processing  "constrained nodes"
    !else

     ! do j = 1, n_vertex_max  
         
       ! node_out(j)=element%vertex(j)

      ! enddo 
   
    ! endif
   
   if (element%n_sons .eq. 0) then 
   
    do i=1,n_vertex_max

      inode         =node_out(i)! element%vertex(i)

      do j=1,n_order+1

        index_ij = (i-1)*(n_order+1) + j     ! index in the ELM matrix

        index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix
 
        mumps_par%rhs(index_large_i) = mumps_par%rhs(index_large_i) + RHS(index_ij)

        do k=1,n_vertex_max

          knode         =node_out(k)! element%vertex(k)

          do l=1,n_order+1

            index_kl = (k-1)*(n_order+1) + l

            index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix

            ilarge = ilarge +1

            mumps_par%irn(ilarge) = index_large_i
            mumps_par%jcn(ilarge) = index_large_k
            mumps_par%A(ilarge)   = ELM(index_ij,index_kl)

          enddo
        enddo
      enddo
    enddo
 endif
  enddo

  nz_AA_old = nz_AA
  nz_AA = ilarge

  zbig = 1.d10

!----------------------- boundary conditions

  do i=1,node_list%n_nodes

    if (node_list%node(i)%boundary .ne. 0) then

      index_i = node_list%node(i)%index(1)  ! base index in the main matrix

      mumps_par%irn(ilarge+1) = index_i
      mumps_par%jcn(ilarge+1) = index_i
      mumps_par%A(ilarge+1)   = zbig
      ilarge = ilarge + 1

      if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then

        index_i = node_list%node(i)%index(2)  ! base index in the main matrix

        mumps_par%irn(ilarge+1) = index_i
        mumps_par%jcn(ilarge+1) = index_i
        mumps_par%A(ilarge+1)   = zbig
        ilarge = ilarge + 1
      endif

      if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3)) then

        index_i = node_list%node(i)%index(3)  ! base index in the main matrix

        mumps_par%irn(ilarge+1) = index_i
        mumps_par%jcn(ilarge+1) = index_i
        mumps_par%A(ilarge+1)   = zbig
        ilarge = ilarge + 1
      endif

    endif
  enddo

  nz_AA_old = nz_AA
  nz_AA     = ilarge

  mumps_par%n  = n_AA
  mumps_par%nz = nz_AA

endif

#ifdef USE_MUMPS
mumps_par%JOB = 6
mumps_par%SYM = 0
mumps_par%icntl(7) = 4

write(*,*) ' mumps : ',mumps_par%n, mumps_par%nz

call DMUMPS(mumps_par)
#else

if (my_id == 0) then
   if (allocated(sparskit_work)) deallocate(sparskit_work)
   allocate(sparskit_work(mumps_par%N + 1))
   print*, "taille du systeme,non zero", n_AA,nz_AA, mumps_par%NZ
   call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)

   nnz = mumps_par%JCN(mumps_par%N+1) - 1
   call pastix_fortran_checkmatrix(check_data, MPI_COMM_SELF, &
        1, 0, 1, mumps_par%N, mumps_par%JCN, mumps_par%IRN, mumps_par%A, -1, 1)
   write (*,*) "nnz", nnz
   mumps_par%NZ = mumps_par%JCN(mumps_par%N+1) - 1
   if (mumps_par%NZ /= nnz ) then
      write (*,*) "associated (mumps_par%IRN)", associated (mumps_par%IRN)
      if (associated (mumps_par%IRN)) deallocate(mumps_par%IRN)
      if (associated (mumps_par%A)  ) deallocate(mumps_par%A)
      allocate(mumps_par%IRN(mumps_par%NZ))
      allocate(mumps_par%A(mumps_par%NZ))
      call pastix_fortran_checkmatrix_end(check_data, &
           1, mumps_par%IRN,mumps_par%A, 1)
   endif

end if

CALL MPI_BCAST(mumps_par%N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(mumps_par%NZ, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
write (*,*) "mumps_par%NZ", mumps_par%NZ
if (my_id /= 0) then
   if (associated(mumps_par%JCN)) deallocate(mumps_par%JCN)
   if (associated(mumps_par%IRN)) deallocate(mumps_par%IRN)
   if (associated(mumps_par%A))   deallocate(mumps_par%A)
   if (associated(mumps_par%rhs)) deallocate(mumps_par%rhs)

   allocate(mumps_par%JCN(mumps_par%N+1))
   allocate(mumps_par%IRN(mumps_par%NZ))
   allocate(mumps_par%A(mumps_par%NZ))
   allocate(mumps_par%RHS(mumps_par%N))
end if

CALL MPI_BCAST(mumps_par%JCN(1), mumps_par%N+1,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(mumps_par%IRN(1), mumps_par%NZ,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(mumps_par%A(1),   mumps_par%NZ,   MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(mumps_par%RHS(1), mumps_par%N,    MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD, ierr)


if (.not. allocated(pastix_perm_vars))  allocate(pastix_perm_vars(mumps_par%n))
if (.not. allocated(pastix_iperm_vars)) allocate(pastix_iperm_vars(mumps_par%n))

pastix_iparm(1)  = 0          ! insert default values
pastix_iparm(2)  = 0          ! initializse
pastix_iparm(3)  = 0

write(*,*) '***********************************'
write(*,*) '* initialise PastiX                *'
write(*,*) '***********************************'
 
pastix_data = 0
 call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
     pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

pastix_iparm(2) = 1
pastix_iparm(3) = 7
pastix_iparm(6) = pastix_iter           ! refinement : max number of iterations

pastix_iparm(7)  = 1                    ! force check

pastix_iparm(31) = pastix_facto
pastix_iparm(35) = pastix_nthrd         ! numthreads : number of threads
pastix_iparm(39) = pastix_rhs           ! right hand side (0 : use RHS)
pastix_iparm(37) = pastix_iluk 
pastix_iparm(41) = pastix_sym

pastix_iparm(42) = pastix_ricar
pastix_iparm(37) = pastix_iluk
pastix_iparm(14) = pastix_amalg

pastix_dparm(6)  = pastix_epsilon    ! error level refinement
pastix_dparm(11) = pastix_pivot      ! pivot threshold?

write(*,*) '***********************************'
write(*,*) '* call PastiX                     *'
write(*,*) '***********************************'

 call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
     pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

#endif
if (my_id .eq. 0) then

  do i=1,node_list%n_nodes
   if(.not. node_list%node(i)%constrained) then 
    index = node_list%node(i)%index(1)
   
!    write(*,'(A,2i5,8e16.8)') 'poisson : ',i,index,mumps_par%RHS(index),mumps_par%RHS(index+1),mumps_par%RHS(index+2),mumps_par%RHS(index+3)

    do k=1,n_order+1

      index = node_list%node(i)%index(k)

      node_list%node(i)%values(i_harm,k,ivar_out) = node_list%node(i)%values(i_harm,k,ivar_out) + (1. - amix) * mumps_par%RHS(index)

    enddo
   endif
  enddo

  !*************************************************************************
  ! Solutions at constrained nodes                                         *
  !*************************************************************************    
! if(refinement==.true.) then  
 		    
  do i = 1, node_list%n_nodes
      
   
       if(node_list%node(i)%constrained) then
	    !write(*,*) ' '
	    !write(*,*) 'Constrained node(poisson)', i
	    !write(*,*)' '
          
	    lambda = node_list%node(i)%ref_lambda
	    mu     = node_list%node(i)%ref_mu
	    index_elm = node_list%node(i)%parent_elem
	    parent(1) = node_list%node(i)%parents(1)
	    parent(2) = node_list%node(i)%parents(2)
	  
	    call basisfunctions(lambda, mu, H, H_s, H_t, H_st)
	    

            do j = 1, n_vertex_max           
                 pr(j) = element_list%element(index_elm)%vertex(j)
           
            end do 
	    
	    dx_ds = 0.
	    dx_dt = 0.
	    dy_ds = 0.
	    dy_dt = 0.
	    d2x_dsdt = 0.
	    d2y_dsdt = 0. 
	          
            Psi = 0.
	    dPsi_ds = 0.
	    dPsi_dt = 0.
	    d2Psi_dsdt = 0.
	    
	      
	         do k = 1, n_vertex_max	        
                    
                      if((pr(k)==parent(1)).or.(pr(k)==parent(2))) then
                          
		      
                           do l = 1, n_order+1

			        dx_ds = dx_ds + node_list%node(pr(k))%x(l,1) * H_s(k,l) 	&
			         * element_list%element(index_elm)%size(k,l)
			        dx_dt = dx_dt + node_list%node(pr(k))%x(l,1) * H_t(k,l) 	&
			         * element_list%element(index_elm)%size(k,l)	  	 

		                dy_ds = dy_ds + node_list%node(pr(k))%x(l,2) * H_s(k,l) 	&
			         * element_list%element(index_elm)%size(k,l)
			        dy_dt = dy_dt + node_list%node(pr(k))%x(l,2) * H_t(k,l) 	&
			         * element_list%element(index_elm)%size(k,l)
				 
			        d2x_dsdt = d2x_dsdt + node_list%node(pr(k))%x(l,1) * H_st(k,l) 	&
			         * element_list%element(index_elm)%size(k,l)	
			        d2y_dsdt = d2y_dsdt + node_list%node(pr(k))%x(l,2) * H_st(k,l) 	&
			         * element_list%element(index_elm)%size(k,l)
		                
				Psi = Psi  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H(k,l)	   &
		                 *element_list%element(index_elm)%size(k,l)				 	 
		                
				dPsi_ds	= dPsi_ds  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_s(k,l)	   &
		                 *element_list%element(index_elm)%size(k,l)
			    
		                dPsi_dt	= dPsi_dt  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_t(k,l)	   &
		                 *element_list%element(index_elm)%size(k,l)
			    
		                d2Psi_dsdt = d2Psi_dsdt  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_st(k,l)	   &
		                 *element_list%element(index_elm)%size(k,l) 
				
				 				  						 	 	  	 				 
                           end do
		
		       end if
			       
	         end do                 
	  
	         h_u = 1
	         h_v = 1
	         h_w = h_u*h_v 
	
            node_list%node(i)%values(i_harm,1,ivar_out) = Psi
            node_list%node(i)%values(i_harm,2,ivar_out)	= (dPsi_ds) /(3.*h_u)
            node_list%node(i)%values(i_harm,3,ivar_out)	= (dPsi_dt) /(3.*h_v)
            node_list%node(i)%values(i_harm,4,ivar_out)	= (d2Psi_dsdt) /(9.*h_w)
         
      end if	
           
  end do

! endif  
  deallocate(mumps_par%irn,mumps_par%jcn,mumps_par%A,mumps_par%rhs)
  !deallocate(pastix_perm_vars,pastix_iperm_vars)
endif


return
end

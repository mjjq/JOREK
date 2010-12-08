subroutine Poisson(my_id,itype,node_list,element_list,bnd_node_list,bnd_elm_list, &
                   ivar_in,ivar_out,i_harm, psi_axis,psi_bnd,xpoint,Z_xpoint,freeboundary,iter)
!---------------------------------------------------------------
! collect the element matrices into one large sparse matrix
! in coordinate format
!---------------------------------------------------------------

use data_structure
use mumps_module

! --- Routine parameters
integer,                  intent(in)    :: my_id
integer,                  intent(in)    :: itype
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list
integer,                  intent(in)    :: ivar_in
integer,                  intent(in)    :: ivar_out
integer,                  intent(in)    :: i_harm
integer,                  intent(in)    :: iter
logical,                  intent(in)    :: xpoint
logical,                  intent(in)    :: freeboundary

! --- local variables
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
type (type_element)      :: element_father
type (type_node)         :: nodes_father(n_vertex_max)
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

real*8   :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1))
real*8   :: zbig, Z_xpoint, psi_axis, psi_bnd, psi_xpoint, R_xpoint, s_xpoint, t_xpoint
real*8   :: R_axis, Z_axis, s_axis, t_axis, amix
integer  :: i_elm_axis, i_elm_xpoint
integer  :: n_AA, nz_AA, nz_AA_old, n_border, ilarge, ife, iv, i,j,k,l
integer  :: n_elements, inode, index_large_i, knode, index_large_k, index_ij, index_kl, index, index_i

real*8, dimension(4,4)	 :: H, H_s, H_t, H_st
real*8			 :: lambda, mu	
real*8			 :: Psi,dPsi_ds,dPsi_dt,d2Psi_dsdt
real*8			 :: dX_ds, dX_dt, dY_ds, dY_dt, d2X_dsdt, d2Y_dsdt, h_u, h_v, h_w
real*8			 :: i_father, inode_father,Index_elm
integer, dimension(n_vertex_max) :: pr
integer, dimension(2)		 :: parent
integer, dimension(n_vertex_max) :: node_out
integer:: nnz, check_data, ierr

write(*,*) '**************************************'
write(*,*) '*            Poisson                 *'
write(*,*) '**************************************'
if (iter .le. 1) then
  write(*,*) ' i_type       : ',itype
  write(*,*) ' n_elements   : ',element_list%n_elements
  write(*,*) ' n_nodes      : ',node_list%n_nodes
  write(*,*) ' freeboundary : ',freeboundary
endif
  
nz_AA = element_list%n_elements * (n_vertex_max * (n_order+1))**2 

n_border = 0
do i=1,node_list%n_nodes
  if (node_list%node(i)%boundary .eq. 1) n_border = n_border+2
  if (node_list%node(i)%boundary .eq. 2) n_border = n_border+2
  if (node_list%node(i)%boundary .eq. 3) n_border = n_border+3
enddo

if ((.not. freeboundary) .or. (itype .ne. -1)) then
  nz_AA = nz_AA + n_border
elseif  (freeboundary .and. (itype .eq. -1)) then
  nz_AA = nz_AA + 128 * bnd_node_list%n_bnd_nodes**2
endif
  
n_AA  = node_list%n_nodes * (n_order+1)

n_AA = 0
do inode = 1, node_list%n_nodes
  n_AA = max(n_AA,node_list%node(inode)%index(4))
enddo

if (iter .le. 1) then
  write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * (n_order+1)
  write(*,*) ' number of boundary nodes: ',n_border
endif
  
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

amix = 0.d0
if (itype .eq. -1) then
  if (freeboundary) then
    amix= 0.95
  else
    amix = 0.5
  endif
endif
    
do ife =1, n_elements

  element = element_list%element(ife)
  i_father= element_list%element(ife)%father

  if (i_father.ne. 0) then
    element_father = element_list%element(i_father)
  endif

  do iv = 1, n_vertex_max

    if (i_father.ne.0) then
      inode_father=element_father%vertex(iv)
      nodes_father(iv) = node_list%node(inode_father)
    endif

    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)

  enddo

  if (itype .eq. -1) then
    
    if (freeboundary) then
      call element_matrix_GS(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
    else
      call element_matrix_GS_perturbation(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
    endif
    
  elseif (itype .eq. -2) then

    call element_matrix_GS_inverse(xpoint,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

  elseif (itype .eq. +2) then

    call element_matrix_Poisson_inverse(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

  else

    call element_matrix_Poisson(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)

  endif

  call Chgmt_node(ife,element,nodes,element_father,nodes_father,ELM,RHS,node_out) ! Processing  "constrained nodes"

  if (element%n_sons .eq. 0) then

    do i=1,n_vertex_max

      inode = node_out(i)! element%vertex(i)

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

  endif      ! element%n_sons

enddo

nz_AA_old = nz_AA
nz_AA = ilarge
mumps_par%nz = nz_AA

zbig = 1.d10

!----------------------- boundary conditions

if (freeboundary .and. (itype .eq. -1)) then

  call vacuum_equil(node_list,bnd_node_list,bnd_elm_list,psi_axis,psi_bnd)
            
else        ! apply fixed boundary conditions

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

mumps_par%n  = n_AA

mumps_par%JOB = 6
mumps_par%SYM = 0
mumps_par%icntl(7) = 4

if (iter .le. 1) write(*,*) ' mumps : ',mumps_par%n, mumps_par%nz

call DMUMPS(mumps_par)

do i=1,node_list%n_nodes

  if(node_list%node(i)%constrained==.false.) then

    index = node_list%node(i)%index(1)

    do k=1,n_order+1

      index = node_list%node(i)%index(k)

!--------------- for equation in perturbation form
      if ((.not. freeboundary) .and. (itype .eq. -1)) then
        node_list%node(i)%deltas(i_harm,k,ivar_out) = mumps_par%RHS(index)
        node_list%node(i)%values(i_harm,k,ivar_out) = node_list%node(i)%values(i_harm,k,ivar_out) &
                                                    + (1.d0 - amix) * mumps_par%RHS(index)
      else
!--------------- for equation on total flux
        node_list%node(i)%deltas(i_harm,k,ivar_out) = node_list%node(i)%values(i_harm,k,ivar_out) - mumps_par%RHS(index)
        node_list%node(i)%values(i_harm,k,ivar_out) = amix * node_list%node(i)%values(i_harm,k,ivar_out) &
                                                    + (1.d0 - amix) * mumps_par%RHS(index)
      endif
    enddo
  endif
enddo

!*************************************************************************
! Solutions at constrained nodes                                         *
!*************************************************************************
! if(refinement==.true.) then

do i = 1, node_list%n_nodes

  if(node_list%node(i)%constrained==.true.) then
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

      if ((pr(k)==parent(1)).or.(pr(k)==parent(2))) then

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
	           * element_list%element(index_elm)%size(k,l)

	  dPsi_ds = dPsi_ds  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_s(k,l)	   &
	          * element_list%element(index_elm)%size(k,l)

          dPsi_dt = dPsi_dt  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_t(k,l)	   &
                  * element_list%element(index_elm)%size(k,l)

          d2Psi_dsdt = d2Psi_dsdt  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_st(k,l)	   &
	          * element_list%element(index_elm)%size(k,l)
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
  
return
end subroutine poisson

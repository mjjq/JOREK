subroutine update_values(my_id,element_list,node_list,RHS)
!-----------------------------------------------------------------------
! subroutine adds the delta_values in RHS to the values in the node_list
!-----------------------------------------------------------------------
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
real*8  :: RHS(*)
real*8, dimension(4,4)	 :: H, H_s, H_t, H_st
real*8, dimension(2,4) 	 :: c, dc_ds, dc_dt, d2c_dsdt					   
real*8			 :: lambda, mu,dif	
real*8, dimension(n_tor) :: Psi, dPsi_ds,dPsi_dt, d2Psi_dsdt
real*8, dimension(n_tor) :: Delt,Delt_ds,Delt_dt,Delt_dsdt
real*8			 :: h_u, h_v, h_w   
integer, dimension(n_vertex_max)  :: Pr
integer, dimension(2)		  :: parent
integer                           :: index_elm,iv,l,i_tor,ivar
integer :: my_id, i, j, k, in, index_node, index

if (my_id .eq. 0) then

  do i = 1, node_list%n_nodes
   if((node_list%node(i)%constrained==.false.) ) then 
    do j=1,n_order+1

      index_node = node_list%node(i)%index(j)

      do k=1,n_var

        do in=1,n_tor

          index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in

          if (index .gt. 0) then

            node_list%node(i)%values(in,j,k) = node_list%node(i)%values(in,j,k) + RHS(index)
            node_list%node(i)%deltas(in,j,k) = RHS(index)

          endif

        enddo

      enddo

    enddo
   
   endif
!   write(*,'(i5,20e12.4)') i,node_list%node(i)%values(1,:,2),node_list%node(i)%values(2,:,2)

  enddo
  !stop
  do i = 1, node_list%n_nodes
   
  if((node_list%node(i)%constrained==.true.) ) then   

            lambda = node_list%node(i)%ref_lambda
	    mu     = node_list%node(i)%ref_mu
	    index_elm = node_list%node(i)%parent_elem
	    parent(1) = node_list%node(i)%parents(1)
	    parent(2) = node_list%node(i)%parents(2)
        
	    call basisfunctions(lambda, mu, H, H_s, H_t, H_st)
	   
            do j = 1, n_vertex_max           
                 Pr(j) = element_list%element(index_elm)%vertex(j)    
            end do 
      
	   
  

	     h_u =1.
	     h_v =1. 
	     h_w =h_u*h_v 

    !***************************************************
    !     update values and deltas                               *
    !***************************************************
   do ivar=1,n_var

            Psi = 0.
	    dPsi_ds = 0.
	    dPsi_dt = 0.
	    d2Psi_dsdt = 0.

            Delt=0.
            Delt_ds=0.
            Delt_dt=0.
            Delt_dsdt=0.  
  
     do i_tor = 1, n_tor     
          do k = 1, n_vertex_max	        
             Pr(k) = element_list%element(index_elm)%vertex(k)    
           if((Pr(k)==parent(1)).or.(Pr(k)==parent(2))) then
            
                 
	 	      
             do l = 1, n_order+1
                         
			       
                        !***************
                        !  Values      *
                        !***************
                 Psi(i_tor) = Psi(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
       
                 dPsi_ds(i_tor) = dPsi_ds(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)

                 dPsi_dt(i_tor) = dPsi_dt(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
           
                 d2Psi_dsdt(i_tor) = d2Psi_dsdt(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar) &
		     * H_st(k,l) *element_list%element(index_elm)%size(k,l)	    
                           
                        !***************
                        !  Deltas      *
                        !***************
 		      
                 Delt(i_tor) = Delt(i_tor) + node_list%node(Pr(k))%deltas(i_tor,l,ivar)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
       
                 Delt_ds(i_tor) = Delt_ds(i_tor) + node_list%node(Pr(k))%deltas(i_tor,l,ivar) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)

                 Delt_dt(i_tor) =  Delt_dt(i_tor) + node_list%node(Pr(k))%deltas(i_tor,l,ivar) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
           
                 Delt_dsdt(i_tor) = Delt_dsdt(i_tor)  + node_list%node(Pr(k))%deltas(i_tor,l,ivar) &
		     * H_st(k,l) *element_list%element(index_elm)%size(k,l)	     
     
				 				 				  						 	 	  	 !print*,"vals",my_id, i,i_tor,Psi(i_tor),dPsi_ds(i_tor),dPsi_dt(i_tor),d2Psi_dsdt(i_tor)
		 !print*,"delt",my_id, i, i_tor,delt(i_tor),delt_ds(i_tor),delt_dt(i_tor),delt_dsdt(i_tor)		 
              end do !(l)
		
            end if
			       
	  end do !(k)                
	
                        !***************
                        !  Values      *
                        !***************
       node_list%node(i)%values(i_tor,1,ivar)	= (Psi(i_tor))
       node_list%node(i)%values(i_tor,2,ivar) 	= (dPsi_ds(i_tor)) / (3.*h_u)
       node_list%node(i)%values(i_tor,3,ivar)	= (dPsi_dt(i_tor)) / (3.*h_v)
       node_list%node(i)%values(i_tor,4,ivar)	= (d2Psi_dsdt(i_tor)) / (9.*h_w)	    
 
      
	                !***************
                        !  Deltas      *
                        !***************
       node_list%node(i)%deltas(i_tor,1,ivar)	= (Delt(i_tor) )
       node_list%node(i)%deltas(i_tor,2,ivar) 	= (Delt_ds(i_tor)) / (3.*h_u)
       node_list%node(i)%deltas(i_tor,3,ivar)	= (Delt_dt(i_tor))/ (3.*h_v)
       node_list%node(i)%deltas(i_tor,4,ivar)	= (Delt_dsdt(i_tor))/ (9.*h_w)
     
    end do!(i_tor)
     
  enddo !(ivar) 


 endif
  
enddo !(i)

endif

call broadcast_nodes(my_id,node_list)

return
end

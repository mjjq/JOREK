!**********************************************************************************
! subroutine refines an existig element if possible otherwise
! returns a list of elements that need to be refined first.
!
! iref : the element to be refined     (input)
! idir : the direction of refinement   (input)
!
! iref_out(1,:) : a list of length n_out with elements that need to be refined before iref
! iref_out(2,:) : the direction for the iref_out list
! n_out    : the number of elements in the lists iref_out, idir_out
! istatus  : an status to indicate the succes or failure
!
! An element can only be refined if its neighbours are of equal size.
! This means a positive value for the neighbour, a negative value indicates
! a larger neighbour.
! A boundary node can always be refined
!
! Refinement can be a split into two (in two possible directions, idir=1,2) 
! or a split into 4. (idir=3)
!
! idir = 1 -> half the element between node (1,2) and (3,4)
! idir = 2 -> half the element between node (2,3) and (4,1)
! idir = 3 -> half the element in both directions
!
!**********************************************************************************


subroutine Refine_Element(node_list, element_list, iref,idir,lambda_ref, mu_ref, iref_out,n_out,istatus,counter)


 use mod_parameters
 use data_structure

 implicit none 

 type (type_node_list)          :: node_list
 type (type_element_list)       :: element_list  

 real*8				:: lambda, mu, lambda_ref, mu_ref

 integer, dimension(2,4)	:: iref_out(2,4)
 integer   			:: iref, iside,istatus,n_out,idir
 integer                        :: ip1,ip2, ip3,ip4,ip5,ip6,ip7,ip8,ip9
 integer                        :: nb1,nb2,inb1,inb2,inb3,inb4,&
                                   ip41,ip34,ip23,ip14,ip12,counter
  


 
  istatus = 999
  if ((iref.lt.0) .or. (iref .gt.element_list%n_elements)) then
       istatus = 1 
       return 
  endif


!-------------- check for sons, fathers can not be refined
  if (element_list%element(iref)%sons(1) .ne. 0) then
       istatus = 3
       return
  endif

!-------------- check valifity of direction of refinement
  if (( idir .lt. 1) .or. (idir .gt. 3)) then
       istatus = 4
       return
  endif
  
  n_out = 0

!-------------- attempt to avoid a deadlock

  if (idir .ne. 3) then
  
       nb1 = idir + 1
       nb2 = mod(nb1 + 1,4) + 1
       inb1 = abs(element_list%element(iref)%neighbours(nb1))

     
       if (inb1 .ne. 0) then
            if (element_list%element(inb1)%sons(1) .ne. 0 ) then
                 idir = 3
            endif
       endif
       inb2 = abs(element_list%element(iref)%neighbours(nb2))
       if (inb2 .ne. 0) then
            if ( element_list%element(inb2)%sons(1) .ne. 0 ) then
                 idir = 3
            endif
       endif
  endif


!-------------- idir = 1 : split between nodes (1,4) and (2,3) corresponding to neighbours 1 and 3 resp.

  call Ref_Check_Neighb_Stat(node_list,element_list,iref,1,iref_out,n_out)
  call Ref_Check_Neighb_Stat(node_list,element_list,iref,3,iref_out,n_out)
  call Ref_Check_Neighb_Stat(node_list,element_list,iref,2,iref_out,n_out)
  call Ref_Check_Neighb_Stat(node_list,element_list,iref,4,iref_out,n_out)
  if (n_out .ne. 0 ) then         ! the element cannot be refined, others needs to be refined first
       istatus = 5
       return                        ! return to do other elements first
  endif

!----------------------------------------------------------------------------------------
! at this point no larger neighbours exist and the element can be refined.
! A neighbour is either of the same size, smaller or does not exist (i.e. a boundary side)
!----------------------------------------------------------------------------------------

  if ((idir .eq. 1) .or. (idir .eq. 3)) then
  
       iside = 1  
       inb1  = element_list%element(iref)%neighbours(iside)  
       ip1   = element_list%element(iref)%vertex(1) 
       ip2   = element_list%element(iref)%vertex(2) 
       lambda =lambda_ref 
       mu     = 0.
 
       
       if ( inb1 .eq. 0 ) then ! if so, this side lies on the boundary  
         
             
            if ( (node_list%node(ip1)%boundary.ne.0) .and. (node_list%node(ip2)%boundary.ne.0) ) then  ! double check if indeed on boundary
                 
		 call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,0,ip12,counter)
		 
            else
                 write(*,*) ' error in data structure :  a zero neighbour is not on the boundary 1 iside',iref,idir,ip1,ip2, iside
            endif
       
       else if ( element_list%element(inb1)%sons(1) .ne. 0 ) then     ! the neighbour has sons meaning the new node already exists
          
            call Ref_Find_Constrained_Node(node_list,element_list,iref,iside,ip12)
            
            if ( node_list%node(ip12)%constrained ) then
                 node_list%node(ip12)%constrained = .false. 
            else
                 write(*,*) ' Problem : this should be a constrained node : ',ip12
            endif
       
       else
           
              call Ref_Add_Node(node_list,element_list,lambda,mu,iref ,iside,+1 ,ip12,counter) ! i12=+1 create a  constrained node
	          
       endif
     
       iside = 3
       inb3  = element_list%element(iref)%neighbours(iside)
      
     
       ip3   = element_list%element(iref)%vertex(3) 
       ip4   = element_list%element(iref)%vertex(4) 
       lambda = lambda_ref
       mu     = 1
         
       if ( inb3 .eq. 0 ) then ! if so, this side lies on the boundary  
  
            if ( node_list%node(ip3)%boundary.ne.0 .and. node_list%node(ip4)%boundary.ne.0) then  ! double check if indeed on boundary
                 call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,0,ip34,counter)
            else
                 write(*,*) ' error in data structure :  a zero neighbour is not on the boundary 2',iref,idir,ip3,ip4
            endif
       
       else if ( element_list%element(inb3)%sons(1) .ne. 0 ) then     ! the neighbour has sons meaning the new node already exists
    
            call Ref_Find_Constrained_Node(node_list,element_list,iref,iside,ip34)    
            if ( node_list%node(ip34)%constrained ) then
                 node_list%node(ip34)%constrained = .false. 
            else
                 write(*,*) ' Problem : this should be a constrained node : ',ip34
            endif
           
	else
    
            call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,+1,ip34,counter)
       
       endif
              
  endif
   
  if ((idir .eq. 2) .or. (idir .eq. 3)) then
  
       iside = 2
       inb2  = element_list%element(iref)%neighbours(iside)
       
       ip2   = element_list%element(iref)%vertex(2) 
       ip3   = element_list%element(iref)%vertex(3) 
       
       lambda = 1.
       mu     = mu_ref
       
       if ( inb2 .eq. 0 ) then ! if so, this side lies on the boundary
       
            if ( node_list%node(ip2)%boundary.ne.0 .and. node_list%node(ip3)%boundary.ne.0 ) then  ! double check if indeed on boundary
                 call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,0,ip23,counter)
              else
                 write(*,*) ' error in data structure :  a zero neighbour is not on the boundary 3 ',iref,idir,ip2,ip3
            endif
	    
         elseif ( element_list%element(inb2)%sons(1) .ne. 0 ) then     ! the neighbour has sons meaning the new node already exists
            call Ref_Find_Constrained_Node(node_list,element_list,iref,iside,ip23)
            if ( node_list%node(ip23)%constrained ) then
                 node_list%node(ip23)%constrained = .false. 
              else
                 write(*,*) ' Problem : this should be a constrained node : ',ip23
            endif
         else
	    lambda = 1.
             mu     = mu_ref
            call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,+1,ip23,counter)
       endif
       
       iside = 4
       inb4  = element_list%element(iref)%neighbours(iside)
       
       
       ip4   = element_list%element(iref)%vertex(4) 
       ip1   = element_list%element(iref)%vertex(1) 
       lambda = 0.
       mu     =mu_ref
      
       if ( inb4 .eq. 0 ) then 
         
            if ( node_list%node(ip4)%boundary.eq.0 .and. node_list%node(ip1)%boundary.eq.0 ) then  ! double check if indeed on boundary "squart"
           
                 call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,-1,ip41,counter)
              else
                 write(*,*) ' error in data structure :  a zero neighbour is not on the boundary 4 ',iref,idir,ip4,ip1
            endif
	    
         else if ( element_list%element(inb4)%sons(1) .ne. 0 ) then     ! the neighbour has sons meaning the new node already exists

            call Ref_Find_Constrained_Node(node_list,element_list,iref,iside,ip41)
            if ( node_list%node(ip41)%constrained ) then
                 node_list%node(ip41)%constrained = .false. 
              else
                 write(*,*) ' Problem : this should be a constrained node : ',ip41
            endif
	    
         else
	 
            call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,+1,ip41,counter)
	    
       endif
    
  endif

  if (idir .eq. 3) then         ! create the new node in the middle
       
       lambda = lambda_ref
       mu = mu_ref
      
       call Ref_Add_Node(node_list,element_list,lambda,mu,iref,iside,-1,ip9,counter)
  endif
 
!------------------------------------------------
! at this point all the new nodes are created 
!------------------------------------------------

  ip5 = ip12
  ip6 = ip34
  ip7 = ip23
  ip8 = ip41
 
  call Ref_Add_Elements(node_list,element_list,iref,idir,lambda_ref,mu_ref, ip5,ip6,ip7,ip8,ip9,istatus)
 
!------------------------------------------------
! finally the neighbours information needs updating
!------------------------------------------------

  call Ref_Update_Neighbours(node_list,element_list,iref,idir)    
 
  istatus = 0
  element_list%element(iref)%nref = abs(idir)
 
  return
  
end subroutine Refine_Element

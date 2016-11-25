!**********************************************************************************
! Subroutine to add the two new elements to the element list
! refined in the direction idir
! P5 and P6 are the indices of the new points of the two elements
! for idir==3 also P7,P8,P9 are required
!**********************************************************************************

subroutine Ref_Add_Elements(node_list,element_list,iref,idir,lambda_ref, mu_ref, &
			    P5,P6,P7,P8,P9,istatus)

 use mod_parameters
 use data_structure

implicit none 

 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list

 real*8, dimension(4,4)	:: H, H_s, H_t, H_st
 real*8, dimension(9)	:: h_u, h_v, h_w
 real*8, dimension(4)	:: signe_hu, signe_hv, coeff_s, coeff_t
 real*8			:: lambda, mu, lambda_ref, mu_ref
 

 integer, dimension(4,4):: local_node
 integer, dimension(4)	:: pr, index_elm 
 integer		:: P1, P2, P3, P4, P5, P6, P7, P8, P9, &
                          new_elm, new_elm1,new_elm2,new_elm3,new_elm4
 integer		:: iref, idir,index_node, istatus
 integer                :: i,j

     !*********************************
     ! Initialisation                 *
     !*********************************
 
  istatus = 999 

  P1 = element_list%element(iref)%vertex(1)  
  P2 = element_list%element(iref)%vertex(2)
  P3 = element_list%element(iref)%vertex(3) 
  P4 = element_list%element(iref)%vertex(4)

  do j = 1, n_vertex_max          
       pr(j) = element_list%element(iref)%vertex(j)       
  end do

  signe_hu(1) =  1.
  signe_hu(2) = -1.
  signe_hu(3) = -1.
  signe_hu(4) =  1.
  
  signe_hv(1) =  1.
  signe_hv(2) =  1.
  signe_hv(3) = -1.
  signe_hv(4) = -1.

     !**************************************************************
     ! create the new son elements and update father/son relations *
     !**************************************************************
  new_elm  = element_list%n_elements
  new_elm1 = new_elm+1
  new_elm2 = new_elm+2
  new_elm3 = new_elm+3
  new_elm4 = new_elm+4
  
  if (new_elm .ge. n_elements_max) then
       write(*,*) 'Add_Elements : AUGMENTER LE NB MAX D''ELEMENTS'
       istatus = 1
       return
  endif

      !****************************
      !  idir= 1                  *
      !****************************
  if (idir .eq. 1) then           
       element_list%element(new_elm1)%vertex(1) = P1;      element_list%element(new_elm1)%vertex(2) = P5
       element_list%element(new_elm1)%vertex(3) = P6;      element_list%element(new_elm1)%vertex(4) = P4
       element_list%element(new_elm1)%father    = iref
       element_list%element(new_elm1)%n_sons    = 0;       element_list%element(new_elm1)%sons   = 0
       element_list%element(new_elm1)%n_gen    = element_list%element(iref)%n_gen + 1

       element_list%element(new_elm2)%vertex(1) = P5;      element_list%element(new_elm2)%vertex(2) = P2  
       element_list%element(new_elm2)%vertex(3) = P3;      element_list%element(new_elm2)%vertex(4) = P6  
       element_list%element(new_elm2)%father    = iref
       element_list%element(new_elm2)%n_sons    = 0;       element_list%element(new_elm2)%sons   = 0 
       element_list%element(new_elm2)%n_gen    = element_list%element(iref)%n_gen + 1
       
       element_list%n_elements = new_elm2
       


       !**************************************************
       ! Bézier coefficients of son elements (idir=1)    *
       !**************************************************
     
       lambda = lambda_ref
       mu = 1.
       
       coeff_s(1) = lambda
       coeff_s(2) = (1.-lambda)

       coeff_t(1) = mu
       coeff_t(2) = mu 
  
     
       !*****************************************************************************
       !Correspondence between nodes of the parent element and nodes of son element *
       !*****************************************************************************
       local_node(1,1) = 1
       local_node(1,2) = 2
       local_node(1,3) = 5
       local_node(1,4) = 4

       local_node(2,1) = 2
       local_node(2,2) = 3
       local_node(2,3) = 6
       local_node(2,4) = 5       
                     
       index_node = 0

       do j = 1, 2
            do i = 1, 3         
                 index_node = index_node + 1      
            end do	    
       end do	    
	    
       index_elm(1) = new_elm1
       index_elm(2) = new_elm2
       	    
       do i = 1, 2
            
            do j = 1, n_vertex_max
                 if( ((i.eq.1).and.(j.eq.2)).or.((i.eq.2).and.(j.eq.1)).or.((i.eq.2).and.(j.eq.4)).or.((i.eq.1).and.(j.eq.3)) )  then 
                     h_u(local_node(i,j))=1.
                     h_v(local_node(i,j))=1.
                     h_w(local_node(i,j))=h_u(local_node(i,j))*h_v(local_node(i,j))
                   else
                     h_u(local_node(i,j))=abs(element_list%element(iref)%size(j,2))
                     h_v(local_node(i,j))=abs(element_list%element(iref)%size(j,3))
                     h_w(local_node(i,j))=abs(element_list%element(iref)%size(j,4))
                   endif
	         element_list%element(index_elm(i))%size(j,1) = 1.
	         element_list%element(index_elm(i))%size(j,2) = coeff_s(i)*signe_hu(j)*h_u(local_node(i,j))
	         element_list%element(index_elm(i))%size(j,3) = coeff_t(i)*signe_hv(j)*h_v(local_node(i,j))
	         element_list%element(index_elm(i))%size(j,4) = coeff_s(i)*signe_hu(j)*signe_hv(j)*coeff_t(i) &
	   					       *h_w(local_node(i,j))
            end do	   			       	  
       end do


       !****************************
       !  idir= 2                 *
       !****************************

    elseif (idir .eq. 2) then


       element_list%element(new_elm1)%vertex(1) = P1;      element_list%element(new_elm1)%vertex(2) = P2
       element_list%element(new_elm1)%vertex(3) = P7;      element_list%element(new_elm1)%vertex(4) = P8
       element_list%element(new_elm1)%father    = iref
       element_list%element(new_elm1)%n_sons    = 0;       element_list%element(new_elm1)%sons   = 0
       element_list%element(new_elm1)%n_gen    = element_list%element(iref)%n_gen + 1
       

       element_list%element(new_elm2)%vertex(1) = P8;      element_list%element(new_elm2)%vertex(2) = P7  
       element_list%element(new_elm2)%vertex(3) = P3;      element_list%element(new_elm2)%vertex(4) = P4  
       element_list%element(new_elm2)%father    = iref
       element_list%element(new_elm2)%n_sons    = 0;       element_list%element(new_elm2)%sons   = 0
       element_list%element(new_elm2)%n_gen    = element_list%element(iref)%n_gen + 1    
          
       element_list%n_elements= new_elm2
       
      !**************************************************
      ! Bézier coefficients of son elements(idir=2)     *
      !**************************************************

       lambda = 1.
       mu = mu_ref
       
       coeff_s(1) = lambda
       coeff_s(2) = lambda

       coeff_t(1) = mu
       coeff_t(2) = 1.- mu 
  
       
       !*****************************************************************************
       !Correspondence between nodes of the parent element and nodes of son element *
       !*****************************************************************************
       local_node(1,1) = 1
       local_node(1,2) = 2
       local_node(1,3) = 4
       local_node(1,4) = 3

       local_node(2,1) = 3
       local_node(2,2) = 4
       local_node(2,3) = 6
       local_node(2,4) = 5       
                     
       index_node = 0

       do j = 1, 3
            do i = 1, 2 
                 index_node = index_node + 1	 
            end do	    
       end do	    
	    
       index_elm(1) = new_elm1
       index_elm(2) = new_elm2
       	    
       do i = 1, 2
            
            do j = 1, n_vertex_max
                if( ((i.eq.1).and.(j.eq.3)).or.((i.eq.1).and.(j.eq.4)).or.((i.eq.2).and.(j.eq.1)).or.((i.eq.2).and.(j.eq.2)) )  then 
                     h_u(local_node(i,j))=1.
                     h_v(local_node(i,j))=1.
                     h_w(local_node(i,j))=h_u(local_node(i,j))*h_v(local_node(i,j))
                   else
                     h_u(local_node(i,j))=abs(element_list%element(iref)%size(j,2))
                     h_v(local_node(i,j))=abs(element_list%element(iref)%size(j,3))
                     h_w(local_node(i,j))=abs(element_list%element(iref)%size(j,4))
                endif

	         element_list%element(index_elm(i))%size(j,1) = 1.
	         element_list%element(index_elm(i))%size(j,2) = coeff_s(i)*signe_hu(j)*h_u(local_node(i,j))
	         element_list%element(index_elm(i))%size(j,3) = coeff_t(i)*signe_hv(j)*h_v(local_node(i,j))
	         element_list%element(index_elm(i))%size(j,4) = coeff_s(i)*signe_hu(j)*signe_hv(j)*coeff_t(i) &
	   					       *h_w(local_node(i,j))

              
            end do	   			       	  
       end do
       
 
     !****************************
     !  idir= 3                  *
     !****************************
    elseif (idir .eq. 3) then

       element_list%element(new_elm1)%vertex(1) = P1;      element_list%element(new_elm1)%vertex(2) = P5
       element_list%element(new_elm1)%vertex(3) = P9;      element_list%element(new_elm1)%vertex(4) = P8

       element_list%element(new_elm1)%father    = iref
       element_list%element(new_elm1)%n_sons    = 0;       element_list%element(new_elm1)%sons   = 0
       element_list%element(new_elm1)%n_gen     = element_list%element(iref)%n_gen + 1
       
 
       element_list%element(new_elm2)%vertex(1) = P5;      element_list%element(new_elm2)%vertex(2) = P2  
       element_list%element(new_elm2)%vertex(3) = P7;      element_list%element(new_elm2)%vertex(4) = P9  

       element_list%element(new_elm2)%father    = iref
       element_list%element(new_elm2)%n_sons    = 0;       element_list%element(new_elm2)%sons   = 0
       element_list%element(new_elm2)%n_gen    = element_list%element(iref)%n_gen + 1
      

       element_list%element(new_elm3)%vertex(1) = P9;      element_list%element(new_elm3)%vertex(2) = P7
       element_list%element(new_elm3)%vertex(3) = P3;      element_list%element(new_elm3)%vertex(4) = P6

       element_list%element(new_elm3)%father    = iref
       element_list%element(new_elm3)%n_sons    = 0;       element_list%element(new_elm3)%sons   = 0
       element_list%element(new_elm3)%n_gen     = element_list%element(iref)%n_gen + 1
       

       element_list%element(new_elm4)%vertex(1) = P8;      element_list%element(new_elm4)%vertex(2) = P9  
       element_list%element(new_elm4)%vertex(3) = P6;      element_list%element(new_elm4)%vertex(4) = P4  

       element_list%element(new_elm4)%father    = iref
       element_list%element(new_elm4)%n_sons    = 0;       element_list%element(new_elm4)%sons   = 0
       element_list%element(new_elm4)%n_gen     = element_list%element(iref)%n_gen + 1           
       element_list%n_elements              = new_elm4



      !**************************************************
      ! Bézier coefficients of son elements  (idir = 3) *
      !**************************************************

       lambda = lambda_ref
       mu = mu_ref
       coeff_s(1) = lambda
       coeff_s(2) = (1.-lambda)
       coeff_s(3) = (1.-lambda)
       coeff_s(4) = lambda
       
       coeff_t(1) = mu
       coeff_t(2) = mu
       coeff_t(3) = (1. - mu)
       coeff_t(4) = (1. - mu)     
       !*****************************************************************************
       !Correspondence between nodes of the parent element and nodes of son element * 
       !*****************************************************************************
       local_node(1,1) = 1
       local_node(1,2) = 2
       local_node(1,3) = 5
       local_node(1,4) = 4 
  
       local_node(2,1) = 2
       local_node(2,2) = 3
       local_node(2,3) = 6
       local_node(2,4) = 5 

       local_node(3,1) = 5
       local_node(3,2) = 6
       local_node(3,3) = 9
       local_node(3,4) = 8 

       local_node(4,1) = 4
       local_node(4,2) = 5
       local_node(4,3) = 8
       local_node(4,4) = 7    
                     
       index_node = 0

       do j = 1, 3
         do i = 1, 3
              index_node = index_node + 1          
         end do	      
       end do	    
	      
       index_elm(1) = new_elm1
       index_elm(2) = new_elm2
       index_elm(3) = new_elm3
       index_elm(4) = new_elm4
       	  

         
        
       do i = 1, 4           
            do j = 1, n_vertex_max
         
                   if(i.ne.j) then 
                     h_u(local_node(i,j))=1.
                     h_v(local_node(i,j))=1.
                     h_w(local_node(i,j))=h_u(local_node(i,j))*h_v(local_node(i,j))
                   else
                     h_u(local_node(i,j))=abs(element_list%element(iref)%size(j,2))
                     h_v(local_node(i,j))=abs(element_list%element(iref)%size(j,3))
                     h_w(local_node(i,j))=abs(element_list%element(iref)%size(j,4))
                   endif
	         element_list%element(index_elm(i))%size(j,1) = 1.
	         element_list%element(index_elm(i))%size(j,2) = coeff_s(i)*signe_hu(j)*h_u(local_node(i,j))
	         element_list%element(index_elm(i))%size(j,3) = coeff_t(i)*signe_hv(j)*h_v(local_node(i,j))
	         element_list%element(index_elm(i))%size(j,4) = coeff_s(i)*signe_hu(j)*signe_hv(j)*coeff_t(i) &
	   					                  *h_w(local_node(i,j))

            end do	   			       	  
       end do
      
  endif
 
       !*************************************
       ! Update  "Father" element           *
       !*************************************

  if ((idir .eq. 1) .or. (idir .eq. 2)) then
       element_list%element(iref)%n_sons    = 2       
       element_list%element(iref)%sons(1)   = new_elm1
       element_list%element(iref)%sons(2)   = new_elm2
       
    elseif (idir .eq. 3) then
    
       element_list%element(iref)%n_sons    = 4
       element_list%element(iref)%sons(1)   = new_elm1
       element_list%element(iref)%sons(2)   = new_elm2
       element_list%element(iref)%sons(3)   = new_elm3
       element_list%element(iref)%sons(4)   = new_elm4
  endif
  
  element_list%element(iref)%nref        = idir  



  istatus = 0

  return

end subroutine Ref_Add_Elements

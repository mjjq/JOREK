!****************************************************************
! Does the refinement of a list of elements
! Input :
!   List_to_be_refined : array with element numbers to be refined
!   n_to_be_refined    : the number of elements in the list
!****************************************************************


subroutine Refine_Elem_List(node_list, element_list,list_to_be_refined,n_to_be_refined)


 use data_structure
 use mod_parameters

 implicit none 

 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list

 real*8	    :: lambda_ref, mu_ref

 integer    ::  list_to_be_refined(n_ref_list), n_to_be_refined, iref_out(2,4)
 integer    ::  iref,idir,n_out,istatus,nref,counter
 integer    ::  i
 logical    ::  critical


  write(*,*) ' '
  write(*,*) '**********************************'
  write(*,*) ' Refinement process	'
  write(*,*) '**********************************'
  write(*,*) ' '

     
  do while ( n_to_be_refined .gt. 0 ) 
  
       if (n_to_be_refined .gt. nref_max ) then
            write(*,*) ' TOO MANY ELEMENTS TO BE REFINED', n_to_be_refined, nref_max 
            return
       endif

       iref = list_to_be_refined(n_to_be_refined)
       !write(*,*) ' iref : ',iref
       counter=0
       if ((iref .le. 0) .or. (iref .gt. element_list%n_elements)) then
            write(*,*) ' Nonsense input in Refine_element_list : ',iref
            n_to_be_refined = n_to_be_refined - 1

         else

            idir = abs(element_list%element(iref)%nref)
            
	    idir = 3
           
	    lambda_ref = 0.5
	    mu_ref = 0.5
            
            call Refine_Element(node_list,element_list,iref,idir,lambda_ref, mu_ref, iref_out,n_out,istatus,counter)
     
            if (istatus .eq. 5 ) then         ! element not refined, others need to be done first
             
                 do i=1,n_out

                      nref = element_list%element(iref_out(1,i))%nref
                      
                      if ( ( nref .eq. 0 ) .or. (abs(nref) .eq. iref_out(2,i)) )  then

                           n_to_be_refined = n_to_be_refined + 1    
                           list_to_be_refined(n_to_be_refined) =  iref_out(1,i)
                           element_list%element(iref_out(1,i))%nref      = -iref_out(2,i)
                         
                        elseif ( ( (nref .eq. -1) .and. (iref_out(2,i) .eq. 2) ) .or. &
                               ( (nref .eq. -2) .and. (iref_out(2,i) .eq. 1) ) ) then 

                           n_to_be_refined = n_to_be_refined + 1    
                           list_to_be_refined(n_to_be_refined) =  iref_out(1,i)
                           element_list%element(iref_out(1,i))%nref      =  -3

                        elseif (nref .eq. -3) then
                           n_to_be_refined = n_to_be_refined + 1    
                           list_to_be_refined(n_to_be_refined) =  iref_out(1,i)
                           element_list%element(iref_out(1,i))%nref      = -iref_out(2,i)
                        else     
                           write(*,*) ' Refine_element_list : One should not normally arrive here',nref,iref_out(2,i)
                           n_to_be_refined = n_to_be_refined - 1  ! this is too forget about the element giving problems
                      endif
                 
                 enddo

              else    ! element has (already) been refined and can be removed from the list
             
                 n_to_be_refined = n_to_be_refined - 1
             
            endif
	 
       endif       
  
  enddo


  return


end subroutine Refine_Elem_List


!*****************************************************************************
! checks whether a neighbour on side iside of element iref is  negative
! if so, this neighbour is added to the list iref_out
! inb_side : the index of the neighbour of iref's neighbour which is iref
!*****************************************************************************

subroutine Ref_Check_Neighb_Stat(node_list, element_list, iref,iside,iref_out,n_out)




 use data_structure

 implicit none 

 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list

 integer, dimension(2,4)	:: iref_out(2,4)
 integer, dimension(1)		:: inb_side_arr
 integer                        :: inb,ifather,inb_abs, inb_side,&
                                   iref,iside,n_out

!**********************************************************************************
! 
!**********************************************************************************

  inb = element_list%element(iref)%neighbours(iside)
   
  if ( inb .lt. 0 ) then       ! the neighbour on this side is larger and needs to be refined first
  
    !--------- if the neighbour of iref is lartger, the neighbours neighbour points to the father of iref!

       ifather      = element_list%element(iref)%father
       inb_abs     = abs(inb)
       inb_side_arr = minloc(element_list%element(inb_abs)%neighbours(1:4), &
                             MASK=(element_list%element(inb_abs)%neighbours(1:4) .eq. ifather) )
       inb_side     = inb_side_arr(1) ! this seems necesary for Compaq F90
      
       if ((inb_side .ge. 1 ) .and. (inb_side .le. 4)) then 

            if (element_list%element(inb_abs)%neighbours(inb_side) .eq. ifather ) then   
                 n_out = n_out + 1
                 iref_out(1,n_out) = abs(inb)
                 iref_out(2,n_out) = mod(inb_side - 1,2) + 1
		
              else
                 write(*,*) ' Error in data structure : Missing neighbour '
                 write(*,*) ' Element ',iref,' has neighbours ',element_list%element(iref)%neighbours
                 write(*,*) ' but neighbour ',inb,' has neighbours ',element_list%element(inb_abs)%neighbours
            endif
	    
         else
	 
            write(*,*) ' Check_Neighbour : Index out of bounds ',inb_side,iref,iside
	    
       endif
   
  endif
  
  return
  
  end subroutine Ref_Check_Neighb_Stat

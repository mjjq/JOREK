!**********************************************************************
! Subroutine to update the neighbours information on the refined 
! edge of the element iref and its sons
!**********************************************************************


subroutine Ref_Update_Neighbours(node_list,element_list,iref,idir)



 use data_structure
 use mod_neighbours

implicit none 

 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list

 integer, dimension(8)	  ::  iside_list, ison_list
 integer                  :: iref,ison,ison1,ison2,ison3,ison4,&
                             nsides,iside,idir, inb,inb_ns, &
                             inb_son,is,ison_side,inb_son_side
 integer                  :: i,j,k


!**********************************************************************
! 
!**********************************************************************

  ison1 = element_list%element(iref)%sons(1)
  ison2 = element_list%element(iref)%sons(2)

  if (idir .eq. 1) then

       element_list%element(ison1)%neighbours(2) = ison2
       element_list%element(ison2)%neighbours(4) = ison1
       element_list%element(ison1)%neighbours(4) = element_list%element(iref)%neighbours(4)
       element_list%element(ison2)%neighbours(2) = element_list%element(iref)%neighbours(2)

    elseif (idir .eq. 2) then

       element_list%element(ison1)%neighbours(3) = ison2
       element_list%element(ison2)%neighbours(1) = ison1
       element_list%element(ison1)%neighbours(1) = element_list%element(iref)%neighbours(1)
       element_list%element(ison2)%neighbours(3) = element_list%element(iref)%neighbours(3)

    elseif (idir .eq. 3) then     ! only the internal boundaries

       ison3 = element_list%element(iref)%sons(3)
       ison4 = element_list%element(iref)%sons(4)
       element_list%element(ison1)%neighbours(2) = ison2
       element_list%element(ison1)%neighbours(3) = ison4
       element_list%element(ison2)%neighbours(3) = ison3
       element_list%element(ison2)%neighbours(4) = ison1
       element_list%element(ison3)%neighbours(1) = ison2
       element_list%element(ison3)%neighbours(4) = ison4
       element_list%element(ison4)%neighbours(2) = ison3
       element_list%element(ison4)%neighbours(1) = ison1

  endif

!------------------------------------- update neighbours neighbours of unrefined side

  if (idir .eq. 1) then 
       iside_list(1:2) = (/ 2, 4 /)
       ison_list(1:2)  = (/ 2, 1 /)
       nsides          = 2
    elseif (idir .eq. 2) then
       iside_list(1:2) = (/ 1, 3 /)
       ison_list(1:2)  = (/ 1, 2 /)
       nsides          = 2
    elseif (idir .eq. 3) then
       nsides          = 0
  endif

  do k=1,nsides
       iside = iside_list(k)
       inb   = element_list%element(iref)%neighbours(iside)
       ison  = element_list%element(iref)%sons(ison_list(k))
       if (inb .gt. 0) then
            where (element_list%element(inb)%neighbours == iref ) element_list%element(inb)%neighbours = ison
            inb_ns = element_list%element(inb)%n_sons
            do i=1,inb_ns
                 inb_son = element_list%element(inb)%sons(i)
                 where (element_list%element(inb_son)%neighbours == -iref ) element_list%element(inb_son)%neighbours = -ison
           enddo
       endif
  enddo

!------------------------------- the refined sides, depends on the refinement of the neighbours

  if (idir .eq. 1) then 
       iside_list(1:4) = (/ 1, 1, 3, 3 /)
       ison_list(1:4)  = (/ 1, 2, 1, 2 /)
       nsides          = 4
    elseif (idir .eq. 2) then
       iside_list(1:4) = (/ 2, 2, 4, 4 /)
       ison_list(1:4)  = (/ 1, 2, 1, 2 /)
       nsides          = 4
    elseif (idir .eq. 3) then
       iside_list(1:8) = (/ 1, 1, 2, 2, 3, 3, 4, 4 /)
       ison_list(1:8)  = (/ 1, 2, 2, 3, 3, 4, 4, 1 /)
       nsides          = 8
  endif

  do k = 1, nsides
  
       iside  = iside_list(k)
       inb    = element_list%element(iref)%neighbours(iside)
       if (inb .gt. 0) then
            inb_ns = element_list%element(inb)%n_sons
            ison   = element_list%element(iref)%sons(ison_list(k))

            if (inb_ns .eq. 0 ) then    ! sons of iref point to larger neighbour
                 element_list%element(ison)%neighbours(iside) = -inb
              else
                 do is = 1,inb_ns
                      inb_son = element_list%element(inb)%sons(is)
                      if ( Neighbours(node_list,element_list%element(ison),element_list%element(inb_son),ison_side,inb_son_side))  then
                           if (ison_side .eq. iside) then
                                element_list%element(ison)%neighbours(ison_side)       = inb_son
                                element_list%element(inb_son)%neighbours(inb_son_side) = ison
                             else
                                write(*,*) ' Update_Neighbours : ,isides do not match'
                           endif
                      endif
                 enddo
            endif
       endif
  enddo

  return

end subroutine Ref_Update_Neighbours

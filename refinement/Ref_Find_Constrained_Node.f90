!*****************************************************************************************
! Subroutine to find the constrained node which is shared by the two
! sons of the iside neighbour of element iref.
! This node has the two nodes from the inb edge of iref as parents
!*****************************************************************************************


subroutine Ref_Find_Constrained_Node(node_list, element_list, iref,iside,inode)

 use mod_parameters
 use data_structure

implicit none
 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list

 integer	:: inb, iv1, iv2, ip1, ip2
 integer	:: iref, iside, inode,n_sons,ison,nv,ip
 integer	:: i, is,iv
 
 

  inb   = element_list%element(iref)%neighbours(iside)

  iv1 = iside                        ! the index of the first vertex of edge iside of element iref
  iv2 = mod(iv1,4)   + 1             ! the second vertex

  ip1 = element_list%element(iref)%vertex(iv1)
  ip2 = element_list%element(iref)%vertex(iv2)

  n_sons = 0
  do i = 1, 4
       if(element_list%element(inb)%sons(i).ne.0) n_sons = n_sons+1
  end do

  do is = 1, n_sons
       ison = element_list%element(inb)%sons(is)
       nv   = n_vertex_max 

       do iv=1,nv
            ip =element_list%element(ison)%vertex(iv)
            if ( ( (node_list%node(ip)%parents(1) .eq. ip1) .and. (node_list%node(ip)%parents(2) .eq. ip2) ) .or. &
               ( (node_list%node(ip)%parents(1) .eq. ip2) .and. (node_list%node(ip)%parents(2) .eq. ip1) ) ) then
                 inode = ip
                 return
            endif
       enddo
  enddo

  write(*,*) ' Find_constrained_node :  NODE NOT fOUND : ',iref,iside

  return
  
end subroutine Ref_Find_Constrained_Node
      

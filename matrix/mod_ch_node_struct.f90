module mod_ch_node_struct
implicit none
contains
subroutine Ch_node_struct(i_element, element,nodes,node_out)
 use mod_parameters
 use data_structure

 type (type_element)   :: element

 type (type_node)      :: nodes(n_vertex_max)

 integer, dimension(n_vertex_max)		::  pos_node_constrained
 integer, dimension(n_vertex_max)               ::  node_out
 integer, dimension(2,2)			:: pos_parent, parent
 
 integer					:: i_element
 integer			 		:: pos1, pos2,n_constrained
 integer                                        :: i, j,k, p					   
 
  do j = 1, n_vertex_max               
       node_out(j)=0 
  enddo 

  n_constrained = 0
  pos_node_constrained = 0
 

  do j = 1, n_vertex_max
      
        if ((.not. nodes(j)%constrained)) then

          node_out(j)=element%vertex(j) 
      
        else
            
	    ! write(*,*) '**************************************'
            !write(*,*) '*   chgt nodes structure: element    *',i_element
            ! write(*,*) '**************************************' 
           
 	  
	    n_constrained = n_constrained + 1
	    pos_node_constrained(n_constrained) = j 
	    
	    parent(n_constrained,1) = nodes(j)%parents(1)
	    parent(n_constrained,2) = nodes(j)%parents(2)
        
	    do k = 1, n_vertex_max   
	         
		 if(element%vertex(k) == parent(n_constrained,1)) then
		      pos_parent(n_constrained,1) = k
		      pos_parent(n_constrained,2) = j		      
		   elseif (element%vertex(k) == parent(n_constrained,2)) then
	              parent(n_constrained,1) = nodes(j)%parents(2)
	              parent(n_constrained,2) = nodes(j)%parents(1)	
		      pos_parent(n_constrained,1) = k
		      pos_parent(n_constrained,2) = j	      
		 endif    
	    enddo    
    
       node_out(j)=parent(n_constrained,2)	    
       endif	    
  enddo
end subroutine Ch_node_struct
end module mod_ch_node_struct

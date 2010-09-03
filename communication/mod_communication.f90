MODULE Communication
  INTERFACE 
     SUBROUTINE vertex_is_local(vertex, islocal)
       INTEGER, INTENT(IN)  :: vertex
       LOGICAL, INTENT(OUT) :: islocal
     END SUBROUTINE vertex_is_local

     SUBROUTINE update_values(my_id,element_list,node_list,RHS)

       type (type_node_list)    :: node_list
       type (type_element_list) :: element_list
       real*8  :: RHS(*)
       integer :: my_id
     END SUBROUTINE update_values

     SUBROUTINE update_deltas(my_id,node_list)
       !---------------------------------------------------------------------
       ! subroutine to create a local list of delta values
       !---------------------------------------------------------------------
       integer               :: my_id
       type (type_node_list) :: node_list
     END SUBROUTINE update_deltas

     SUBROUTINE split_brodcast(type,MPI_COMM_N)
       !
       ! Split MPI_BCAST if MPI buffer beyond 2Go
       ! 
       !  
       CHARACTER*8           :: type
       INTEGER               :: MPI_COMM_N
     END SUBROUTINE split_brodcast

  END INTERFACE
END MODULE Communication

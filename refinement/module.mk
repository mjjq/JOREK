DIR = refinement

ALL_BINARIES_SRC := $(ALL_BINARIES_SRC)			\
	$(DIR)/neighbours.f90

JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 		\
	$(DIR)/Ref_Add_Elements.f90		\
	$(DIR)/Ref_Add_Node.f90			\
	$(DIR)/Ref_boundary_node.f90		\
	$(DIR)/Ref_Check_Neighb_Stat.f90	\
	$(DIR)/Ref_Find_Constrained_Node.f90	\
	$(DIR)/Refine_Element.f90		\
	$(DIR)/Refine_Elem_List.f90		\
	$(DIR)/Ref_Update_Neighbours.f90	\
	$(DIR)/Ref_Update_Index.f90		\

DIR=refinement
JOREK2_MAIN_SRC := $(JOREK2_MAIN_SRC) 		\
	$(DIR)/neighbours.f90			\
	$(DIR)/Ref_Add_Elements.f90		\
	$(DIR)/Ref_Add_Node.f90			\
	$(DIR)/Ref_boundary_node.f90		\
	$(DIR)/Ref_Check_Neighb_Stat.f90	\
	$(DIR)/Ref_Find_Constrained_Node.f90	\
	$(DIR)/Refine_Element.f90		\
	$(DIR)/Refine_Elem_List.f90		\
	$(DIR)/Ref_Update_Neighbours.f90	\
	$(DIR)/Ref_Update_Index.f90		\

JOREK2_POINCARE_SRC := $(JOREK2_POINCARE_SRC)	\
	$(DIR)/neighbours.f90

JOREK2_CONNECTION2_SRC := $(JOREK2_CONNECTION2_SRC)	\
	$(DIR)/neighbours.f90

JOREK2_STRIKES_SRC := $(JOREK2_STRIKES_SRC)	\
	$(DIR)/neighbours.f90

JOREK2VTK_SRC := $(JOREK2VTK_SRC)

JOREK2FLVTK_SRC := $(JOREK2FLVTK_SRC)		\
	$(DIR)/neighbours.f90

JOREK2VTK3D_SRC := $(JOREK2VTK3D_SRC)

JOREK2_DIAGNO_SRC := $(JOREK2_DIAGNO_SRC)

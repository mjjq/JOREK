module refinement_module
contains 
! #include rather than include for easiest dependency computation
#include "Ref_Add_Elements.f90"
#include "Ref_Add_Node.f90"
#include "Ref_boundary_node.f90"
#include "Ref_Check_Neighb_Stat.f90"
#include "Ref_Find_Constrained_Node.f90"
#include "Refine_Element.f90"
#include "Refine_Elem_List.f90"
#include "Ref_Update_Neighbours.f90"
#include "Ref_Update_Index.f90"
end module refinement_module

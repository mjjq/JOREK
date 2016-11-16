module mod_elt_matrix_fft
contains
subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
implicit none
 
type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), pointer  :: ELM
real*8, dimension (:)  , pointer  :: RHS
integer, intent(in) :: tid, xcase2

real*8     :: minRad, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
logical    :: xpoint2

real*8, dimension(:,:,:) , pointer :: ELM_p
real*8, dimension(:,:,:) , pointer :: ELM_n
real*8, dimension(:,:,:) , pointer :: ELM_k
real*8, dimension(:,:,:) , pointer :: ELM_kn
real*8, dimension(:,:)   , pointer :: RHS_p
real*8, dimension(:,:)   , pointer :: RHS_k 

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_p
real*8, dimension(:,:,:,:) , pointer :: eq_ss, eq_st, eq_tt   
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t

eq_g    => thread_struct(tid)%eq_g   
eq_s    => thread_struct(tid)%eq_s   
eq_t    => thread_struct(tid)%eq_t   
eq_p    => thread_struct(tid)%eq_p   
eq_ss   => thread_struct(tid)%eq_ss  
eq_st   => thread_struct(tid)%eq_st  
eq_tt   => thread_struct(tid)%eq_tt  
delta_g => thread_struct(tid)%delta_g
delta_s => thread_struct(tid)%delta_s
delta_t => thread_struct(tid)%delta_t

ELM_p  => thread_struct(tid)%ELM_p 
ELM_n  => thread_struct(tid)%ELM_n 
ELM_k  => thread_struct(tid)%ELM_k
ELM_kn => thread_struct(tid)%ELM_kn
RHS_p  => thread_struct(tid)%RHS_p 
RHS_k  => thread_struct(tid)%RHS_k 
return
end subroutine element_matrix_fft
end module mod_elt_matrix_fft
      

!> Contains functions to calculate poloidal currents and the associated FF' 
!!
!!  * Poloidal currents are calculated from teh equilibrium assumption JxB=\grad p 
module mod_poloidal_currents
  
  use constants
  use mod_parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use tr_module
  use phys_module
  use mod_interp
  
  implicit none
  
  private
  
  public :: J_pol 
  
  contains
  
  !!-------------------------------------------------------------------
  !> Calculates poloidal currents from equilibrium assumption
  !!
  !!                  JxB = \grad p
  !! 
  !!    which gives in the JOREK coordinate system and variables
  !!
  !!      J_R = (-j * B_R - R * dp/dZ) / F0
  !!      J_Z = (-j * B_Z + R * dp/dR) / F0
  !!-------------------------------------------------------------------
  subroutine J_pol(node_list, element_list, i_elm, s, t, i_plane, JR, JZ, axisym)

    implicit none

    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    
    integer, intent(in)       :: i_elm, i_plane ! element index / toroidal plane 
    real*8,  intent(in)       ::   s,  t        ! s-t local coordinates
    real*8,  intent(inout)    ::  JR, JZ        ! output current density
    logical, intent(in)       :: axisym         ! if true, only calculated axisymmetric component

    ! --- local variables    
    real*8     :: Psi,Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt
    real*8     :: ZJ ,ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt
    real*8     :: RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
    real*8     :: Ti0,Ti0_s,Ti0_t,Ti0_st,Ti0_ss,Ti0_tt
    real*8     :: Te0,Te0_s,Te0_t,Te0_st,Te0_ss,Te0_tt
    real*8     :: T0,T0_s,T0_t,T0_st,T0_ss,T0_tt
    real*8     :: R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
    real*8     :: xjac, psi_x, psi_y, T0_x, T0_y, RHO_x, RHO_y
    real*8     :: BR, BZ, P0_R, P0_Z, zj_sum, rho_sum, T0_sum
    integer    :: i_tor
   
    call interp_RZ(node_list,element_list,i_elm,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
    xjac  = R_s * Z_t - R_t * Z_s

    psi_x  = 0.d0;   psi_y = 0.d0;    T0_x = 0.d0;     T0_y = 0.d0 
    RHO_x  = 0.d0;   RHO_y = 0.d0;  zj_sum = 0.d0;  rho_sum = 0.d0
    T0_sum = 0.d0 

    do i_tor=1, n_tor

      if ( ( i_tor > 1 ) .and. axisym  ) cycle ! Just include the n=0 mode

      call interp(node_list,element_list,i_elm,1,i_tor,s,t,Psi,Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt)
      call interp(node_list,element_list,i_elm,3,i_tor,s,t,ZJ ,ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt)
      call interp(node_list,element_list,i_elm,5,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
      if (jorek_model .eq. 400) then
        call interp(node_list,element_list,i_elm,6,1,s,t,Ti0,Ti0_s,Ti0_t,Ti0_st,Ti0_ss,Ti0_tt)
        call interp(node_list,element_list,i_elm,8,1,s,t,Te0,Te0_s,Te0_t,Te0_st,Te0_ss,Te0_tt)
        T0    = Ti0   + Te0
        T0_s  = Ti0_s + Te0_s
        T0_t  = Ti0_t + Te0_t
      else
        call interp(node_list,element_list,i_elm,6,i_tor,s,t,T0,T0_s,T0_t,T0_st,T0_ss,T0_tt) 
      endif

      zj_sum  = zj_sum   +  ZJ * HZ(i_tor,i_plane)
      rho_sum = rho_sum  + RHO * HZ(i_tor,i_plane)
      T0_sum  = T0_sum   +  T0 * HZ(i_tor,i_plane)

      if (abs(xjac) > 1.d-6) then 

        psi_x = psi_x + (   Z_t * Ps_s - Z_s * Ps_t )   / xjac * HZ(i_tor,i_plane)
        psi_y = psi_y + ( - R_t * Ps_s + R_s * Ps_t )   / xjac * HZ(i_tor,i_plane)

        T0_x  = T0_x  + (   Z_t * T0_s - Z_s * T0_t )   / xjac * HZ(i_tor,i_plane)
        T0_y  = T0_y  + ( - R_t * T0_s + R_s * T0_t )   / xjac * HZ(i_tor,i_plane)
 
        RHO_x = RHO_x + (   Z_t * RHO_s - Z_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
        RHO_y = RHO_y + ( - R_t * RHO_s + R_s * RHO_t ) / xjac * HZ(i_tor,i_plane)       

      endif
       
    enddo   ! n_tor loop     

    P0_R = T0_sum * RHO_x + T0_x * RHO_sum
    P0_Z = T0_sum * RHO_y + T0_y * RHO_sum
    BR   =  psi_y/R
    BZ   = -psi_x/R

    JR   = ( -ZJ_sum * BR - R * P0_Z ) / F0
    JZ   = ( -ZJ_sum * BZ + R * P0_R ) / F0

  end subroutine J_pol 
  
end module mod_poloidal_currents

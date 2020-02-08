subroutine find_crossing(node_list,element_list,surface_list,j_surf,R_c,Z_c, &
                         R_out,Z_out,ielm_flux,r_flux,s_flux,t_tht,ifail, gofast)

  use data_structure
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  type (type_surface_list), intent(in)          :: surface_list
  integer,                  intent(in)          :: j_surf
  real*8,                   intent(inout)       :: R_c(4),   Z_c(4)
  real*8,                   intent(inout)       :: R_out,    Z_out
  real*8,                   intent(inout)       :: r_flux,s_flux,t_tht
  integer,                  intent(inout)       :: ifail
  logical,                  intent(in)          :: gofast

  ! --- Local variables
  integer :: k
  integer :: ielm_flux
  
  if ((R_c(1) .eq. R_c(3)) .and. (Z_c(1) .eq. Z_c(3))) then
    ifail = 9
    return
  endif
  
  ielm_flux = 0
  
  do k=1,surface_list%flux_surfaces(j_surf)%n_pieces
    
    call find_crossing_on_surface_piece(node_list,element_list,surface_list%flux_surfaces(j_surf), &
                                        k, R_c,Z_c, R_out,Z_out, r_flux,s_flux, t_tht, ifail, gofast)
    if (ifail .eq. 0) then
      ielm_flux = surface_list%flux_surfaces(j_surf)%elm(k)
      return
    endif
  
  enddo
  
  if (ielm_flux .eq. 0) ifail = 99
  
  !write(*,'(A,8e16.8)') ' crossing wrong exit ',x,errx,errf
  
  return

end

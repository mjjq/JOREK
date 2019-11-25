!> Finds the value of psi which contour encloses a given area
!!
!! It's used to define psi_bnd during the freeboundary equilibrium
!! for limiter cases. Psi_bnd is therefore defined in each iteration  
!! as the psi which countour encloses the area of the target equilibrium.
!! The convergence of this method is very slow but allows to run the
!! equilibrium without feedback on the radial position.

module mod_iterate2area

  use data_structure
  
  implicit none
  
  contains
  
  
  !------------------------------------------------------------------
  !> Iterative algorithm to search for psi_bnd with the target area
  !------------------------------------------------------------------
  subroutine iterate2area(node_list,element_list, psi_axis, psi_lim, xpoint2, xcase2, area_ref, psi_bnd)
    
    use phys_module, only : R_limiter, Z_limiter
        
    implicit none
    
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    
    integer,  intent(in)     :: xcase2
    real*8,   intent(in)     :: psi_axis, psi_lim, area_ref
    real*8,   intent(inout)  :: psi_bnd
    logical,  intent(in)     :: xpoint2
    
    ! --- Local variables
    real*8                   :: psi_in, psi_med, psi_out, area, tolerance, Rlim, Zlim
    integer                  :: i, n_max_iter
    
    n_max_iter = 60
    tolerance  = 1.d-6
    
    ! --- Define search limits for psi
    psi_in  = psi_axis     
    psi_out = psi_lim
    
    ! --- Search target area between the given limits
    do i=1, n_max_iter
         
      psi_med = (psi_out + psi_in) * 0.5d0            
      call area_inside_flux_contour(node_list,element_list, xpoint2, xcase2, psi_med, area, Rlim, Zlim)
      
      if (abs(area) > abs(area_ref) ) then                                                             
        psi_out = psi_med                                                                  
      else                                                                          
        psi_in  = psi_med                                                                   
      endif
      
      if(abs((area-area_ref)/area_ref) < tolerance) then
        write(*,*) 'Psi_bnd found for target area = ', psi_med
        exit
      endif
            
    enddo
    
    psi_bnd = psi_med
    R_limiter(1) = Rlim
    Z_limiter(1) = Zlim
    write(*,'(A,2es14.6)') ' area_found, area_ref ', area, area_ref
  
  end subroutine iterate2area
  
       
       
       
       
  
  !------------------------------------------------------------------
  !> Calculates the poloidal cross section area inside a flux contour
  !------------------------------------------------------------------
  subroutine area_inside_flux_contour(node_list,element_list, xpoint2, xcase2, psi, area, Rlim, Zlim)
    
    use tr_module
    use mod_interp, only: interp_RZ
    
    implicit none
    
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    
    integer,  intent(in)     :: xcase2
    real*8,   intent(in)     :: psi
    real*8,   intent(inout)  :: area, Rlim, Zlim
    logical,  intent(in)     :: xpoint2
    
    ! --- local variables    
    integer             :: i, j, n_bnd_points, n_points_piece, n_pieces, min_pos
    integer             :: i_elm, ip, k, node1, node2, node3, node4, counter
    real*8              :: t, rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi
    real*8              :: R0, Z0
    real*8, allocatable :: Rbnd(:), Zbnd(:), Rnew(:), Znew(:), angle(:)
    type (type_surface_list) :: surface_list
    
    n_points_piece = 40
    
    ! --- Define and find flux surface
    surface_list%n_psi =1
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    surface_list%psi_values(1) = psi
    call find_flux_surfaces(999,xpoint2,xcase2,node_list,element_list,surface_list)  
    
    n_pieces = surface_list%flux_surfaces(1)%n_pieces
    
    if (n_pieces < 2) return  ! stops if find_flux_surfaces fails
    n_bnd_points = n_pieces * n_points_piece
    
    if (allocated(Rbnd))  deallocate(Rbnd)
    if (allocated(Zbnd))  deallocate(Zbnd)
    if (allocated(Rnew))  deallocate(Rnew)
    if (allocated(Znew))  deallocate(Znew)
    if (allocated(angle)) deallocate(angle)
    
    allocate(Rbnd(n_bnd_points),Zbnd(n_bnd_points),Rnew(n_bnd_points),Znew(n_bnd_points),angle(n_bnd_points))
    
    ! --- Find n_bnd_points along the flux surface, Rbnd, Zbnd
    
    do k=1, n_pieces
    
        i_elm = surface_list%flux_surfaces(1)%elm(k)
    
        node1 = element_list%element(i_elm)%vertex(1)
        node2 = element_list%element(i_elm)%vertex(2)
        node3 = element_list%element(i_elm)%vertex(3)
        node4 = element_list%element(i_elm)%vertex(4)
    
        rr1  = surface_list%flux_surfaces(1)%s(1,k)
        drr1 = surface_list%flux_surfaces(1)%s(2,k)
        rr2  = surface_list%flux_surfaces(1)%s(3,k)
        drr2 = surface_list%flux_surfaces(1)%s(4,k)
    
        ss1  = surface_list%flux_surfaces(1)%t(1,k)
        dss1 = surface_list%flux_surfaces(1)%t(2,k)
        ss2  = surface_list%flux_surfaces(1)%t(3,k)
        dss2 = surface_list%flux_surfaces(1)%t(4,k)
    
        do ip = 1, n_points_piece
    
          t = -1. + 2.*float(ip-1)/float(n_points_piece-1)
          counter = ip + (k-1)*n_points_piece
    
          call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
          call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)
    
          call interp_RZ(node_list,element_list,i_elm,ri,si,Rbnd(counter),Zbnd(counter))
        enddo
    enddo
    
    if (allocated(surface_list%psi_values))        call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    if (allocated(surface_list%flux_surfaces))     deallocate(surface_list%flux_surfaces)
    
    ! --- Sort points by polar angle
    R0 = (maxval(Rbnd) + minval(Rbnd) ) * 0.5d0
    Z0 = (maxval(Zbnd) + minval(Zbnd) ) * 0.5d0
    do i = 1, n_bnd_points
      angle(i) = atan2(Zbnd(i)-Z0, Rbnd(i)-R0)
    enddo
    
    do i = 1, n_bnd_points
      min_pos = minloc(angle,1)
      Rnew(i) = Rbnd(min_pos)
      Znew(i) = Zbnd(min_pos)
      angle(min_pos) = 100.d0
    enddo
    
    Rnew(n_bnd_points) = Rnew(1) ! First and last point should be identical (for line integral)
    Znew(n_bnd_points) = Znew(1)
    
    !--- Calculate area inside the closed curve (using Green's theorem)
    area = 0.d0
    do i = 1, n_bnd_points-1
      area =  area + (Rnew(i)+Rnew(i+1))*(Znew(i+1)-Znew(i)) - (Znew(i)+Znew(i+1))*(Rnew(i+1)-Rnew(i))
    enddo
    area = 0.25d0*area    
    
    ! --- Save limiter points
    Rlim = Rnew(1)
    Zlim = Znew(1)
      
    deallocate(Rbnd, Zbnd, Rnew, Znew, angle)    

  end subroutine area_inside_flux_contour
  
    
end module mod_iterate2area

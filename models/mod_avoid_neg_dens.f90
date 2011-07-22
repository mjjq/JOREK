!> Increase particle diffusivity around positions where negative densities would develop otherwise.
!! 
!! Negative particle densities can occur in JOREK in the vicinity of large gradients as they occur
!! around ballooning-fingers of ELMs as the numerical scheme is not monotonicity-preserving. This
!! module provides not a clean numerical solution to this problem but a "quick and dirty" solution
!! which, however, already has shown to be useful for some cases.
!!
!! author: Matthias Hoelzl, IPP Garching, 2011.
module avoid_neg_dens
  
  use parameters, only: n_plane
  
  implicit none
  
  private
  public identify_dens_problems, diffusivity_factor, output_problem_pos
  
  !> Use which method for weighting problem positions by rho/rho1 (1: none, 2: linear, 3: quadratic, 4: lin+quadratic)
  integer, parameter :: weighting_type = 3
  
  !> Use which function for increasing the diffusivity (1: Gauss, 2: Tanh, 3: Double-Gauss)
  integer, parameter :: function_type = 3
  
  !> Use which combining method (1: Max, 2: sqrt(sum(a_i^2)))
  integer, parameter :: combining_type = 2
  
  !> Factor by which the diffusivity is locally increased.
  real*8,  parameter :: factor = 120.d0
  
  !> Width around a position with density below vacuum density with increased diffusivity.
  real*8,  parameter :: width  = 0.02d0
  
  !> Problem positions are detected in the region psi_n < psi_n_limit only.
  real*8,  parameter :: psi_n_limit = 1.06d0
  
  !> Maximum number of positions with density below vacuum density per toroidal plane
  integer, parameter :: MAX_PROBLEM_POS = 20001
  
  !> Data structure for storing a position with density below vacuum density: rho < rho_1
  type :: type_problem_pos
    real*8 :: R             !< R-position
    real*8 :: Z             !< Z-position
    real*8 :: rho_over_rho1 !< Ratio between local density, rho, and vacuum density, rho_1
  end type type_problem_pos
  
  !> List of positions with density below vacuum density
  type(type_problem_pos) :: problem_pos(n_plane, MAX_PROBLEM_POS)
  
  !> Number of positions with density below vacuum density
  integer                :: num_problem_pos(n_plane)
  
  
  
  contains
  
  
  
  !> Identify all positions with rho < rho_1
  subroutine identify_dens_problems()
    
    use parameters,        only: n_vertex_max, n_order, n_tor
    use data_structure,    only: type_element, type_node
    use nodes_elements,    only: node_list, element_list
    use gauss,             only: n_gauss
    use basis_at_gaussian, only: H, HZ
    use phys_module,       only: rho_1
    
    type(type_element) :: element
    type(type_node)    :: node
    integer :: mp, ielm, ms, mt, i, j, in, iv
    integer, parameter :: k_psi = 1, k_rho = 5
    real*8  :: R, Z, psi, rho, psi_n
    real*8  :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint
    real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
    integer :: i_elm_xpoint, i_elm_axis, ifail
    
    num_problem_pos(:) = 0
    
    call find_xpoint(0,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,ifail)
    call find_axis(0,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
    
    do mp = 1, n_plane
      do ielm = 1, element_list%n_elements
        element = element_list%element(ielm)
        do ms = 1, n_gauss
          do mt = 1, n_gauss
            R   = 0.d0
            Z   = 0.d0
            psi = 0.d0
            rho = 0.d0
            
            do i = 1, n_vertex_max
              iv   = element%vertex(i)
              node = node_list%node(iv)
              do j = 1, n_order + 1
                
                R = R + node%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
                Z = Z + node%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
                
                do in = 1, n_tor
                  
                  psi = psi + node%values(in,j,k_psi) * element%size(i,j) * H(i,j,ms,mt) * HZ(in,mp)
                  rho = rho + node%values(in,j,k_rho) * element%size(i,j) * H(i,j,ms,mt) * HZ(in,mp)
                  
                end do
                
              end do
            end do
            
            psi_n = ( psi - psi_axis ) / ( psi_xpoint - psi_axis )
            
            if ( ( rho < 0.999d0 * rho_1 ) .and. ( Z > Z_xpoint ) .and. ( psi_n < psi_n_limit ) ) then
              num_problem_pos(mp) = num_problem_pos(mp) + 1
              problem_pos(mp,num_problem_pos(mp))%R = R
              problem_pos(mp,num_problem_pos(mp))%Z = Z
              problem_pos(mp,num_problem_pos(mp))%rho_over_rho1 = rho / rho_1
            end if
            
          end do
        end do
      end do
    end do
    
  end subroutine identify_dens_problems
  
  
  
  !> Factor by which the diffusivity is increased locally at the given position.
  real*8 recursive function diffusivity_factor(R, Z, iplane)
    real*8,  intent(in) :: R, Z
    integer, intent(in) :: iplane
    
    integer :: ipos
    real*8  :: distance, rho_over_rho1, weighting, function_value
    
    diffusivity_factor = 0.d0
    
    do ipos = 1, num_problem_pos(iplane)
    
      distance      = sqrt( (R-problem_pos(iplane,ipos)%R)**2 + (Z-problem_pos(iplane,ipos)%Z)**2 )
      rho_over_rho1 = problem_pos(iplane,ipos)%rho_over_rho1
      
      if ( weighting_type == 1 ) then
        weighting = 1.d0
      else if ( weighting_type == 2 ) then
        weighting = ( 1.d0 - rho_over_rho1 )
      else if ( weighting_type == 3 ) then
        weighting = ( 1.d0 - rho_over_rho1 )**2
      else if ( weighting_type == 4 ) then
        weighting = ( 1.d0 - rho_over_rho1 ) + ( 1.d0 - rho_over_rho1 )**2
      end if
      
      if ( function_type == 1 ) then
        function_value = exp( - distance**2 / width**2 )
      else if ( function_type == 2 ) then
        function_value = 0.5d0 * ( 1.d0 + Tanh( (distance-width) / (-0.4*width) ) )
      else if ( function_type == 3 ) then
        function_value = 0.7d0 * exp( - distance**2 / width**2 ) + 0.3d0 * exp( - distance**2 / (2.d0*width)**2 )
      end if
      
      if ( combining_type == 1 ) then
        diffusivity_factor = max( diffusivity_factor, (factor-1.d0) * function_value * weighting )
      else if ( combining_type == 2 ) then
        diffusivity_factor = diffusivity_factor + ( (factor-1.d0) * function_value * weighting )**2
      end if
      
      
    end do
    
    if ( combining_type == 1 ) then
      diffusivity_factor = diffusivity_factor + 1.d0
    else if ( combining_type == 2 ) then
      diffusivity_factor = sqrt(diffusivity_factor) + 1.d0
    end if
    
  end function diffusivity_factor
  
  
  
  !> Output all positions with rho < rho_1 to an ascii file problem_pos.dat.
  subroutine output_problem_pos()
    
    integer :: iplane, ipos
    real*8  :: min_rho_over_rho1
    
    min_rho_over_rho1 = 1.d0
    
    write(*,'(a,i6)') 'Number of positions with rho < rho_1:', sum(num_problem_pos(1:n_plane))
    
    open(42, file='problem_pos.dat', status='replace', action='write')
    do iplane = 1, n_plane
      if ( num_problem_pos(iplane) == 0 ) write(42,*) '1.6 NaN NaN'
      do ipos = 1, num_problem_pos(iplane)
        write(42,'(3es15.7)') problem_pos(iplane,ipos)%R, problem_pos(iplane,ipos)%Z, problem_pos(iplane,ipos)%rho_over_rho1
        min_rho_over_rho1 = min( min_rho_over_rho1, problem_pos(iplane,ipos)%rho_over_rho1 )
      end do
      write(42,*)
      write(42,*)
    end do
    close(42)
    
    write(*,'(a,es15.7)') 'Worst rho_over_rho1: ', min_rho_over_rho1
    
  end subroutine output_problem_pos
  
  
  
end module avoid_neg_dens

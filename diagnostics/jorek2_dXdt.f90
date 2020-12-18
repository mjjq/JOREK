program jorek2_dXdt

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp
use pellet_module
use mod_bootstrap_functions
use corr_neg
use mod_import_restart
use equil_info
use mod_boundary
use mod_interp

implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list 

allocate(node_list)
allocate(element_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                              ! define the basis functions at the Gaussian points

call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

minRad = 0.0
if (bootstrap) then
  call bootstrap_find_minRad(node_list, element_list, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd)
  call bootstrap_get_q_and_ft_splines(node_list, element_list, ES%psi_axis, ES%psi_xpoint, ES%R_xpoint, ES%Z_xpoint)
  call bootstrap_get_averaged_j_spline(node_list, element_list, ES%psi_axis, ES%psi_xpoint, ES%R_xpoint, ES%Z_xpoint)
endif

do i=1,element_list%n_elements
  
  a(:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss
  
          x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(j,1) * element%size(i,j) * H_t(i,j,ms,mt)
  
          x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(j,1) * element%size(i,j) * H_tt(i,j,ms,mt)
  
          y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H(i,j,ms,mt)
          y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(j,2) * element%size(i,j) * H_t(i,j,ms,mt)
  
          y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
          y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_st(i,j,ms,mt)
          y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
  
        end do
      end do
  
      do ms=1, n_gauss
        do mt=1, n_gauss
          do k=1,n_var
            do in=1,n_tor
              do mp=1,n_plane
                eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
                eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
                eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
  
                eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
                eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
                eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)
  
                delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
                delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
                delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)
              enddo
  
            enddo
  
          enddo
  
        enddo
      enddo
    enddo
  enddo

  do ms = 1, n_gauss
    do mt = 1, n_gauss
      do XXX = 1, n_plane
        ! diff, sources, ...
        xjac  = R_s * Z_t - R_t * Z_s
          do mtor = 1, n_tor
            do i = ! vertices
              do j = ! dofs
  
                v   =  H(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
                v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HHZ(im,mp)

                a(0) = a(0) + v * BigR * xjac * tstep ! needed as \int dV (normalization w.r.t. the volume)
                a(1) = a(1) + v * BigR * (particle_source(ms,mt) + source_pellet)                      * xjac * tstep
              
  end do
  !write out (incl coord and Psi_n)
end do

end program jorek2_dXdt

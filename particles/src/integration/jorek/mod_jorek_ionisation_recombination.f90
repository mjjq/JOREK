module mod_jorek_ionisation_recombination
use mod_ionisation_recombination
implicit none
private
public particle_ion_rec

contains
!> Calculate the new charge of a particle after a specific time
subroutine update_particle_charge(this, sim, particle)
use mod_particle_types
use data_structure
use mod_constants
use phys_module
implicit none
interface
  pure subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
    import :: type_node_list, type_element_list
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ
end interface

type(type_node_list), intent(in)   :: node_list
type(type_element_list), intent(in):: element_list
type(type_particle), intent(inout) :: particle
type(type_ADF11_all), intent(in)   :: ad !< ADF11 datatype
real*8, intent(in)                 :: timestep !< Timestep, in s

real*8, dimension(2) :: P, P_s, P_t, P_phi
real*8               :: R, R_s, R_t, Z, Z_s, Z_t
integer :: i_var(2)
real*8  :: temperature, n_e, T_e

! Get the indices of temperature and density
#ifdef FULLMHD
i_var = (/8,7/) ! T, rho
#else
#ifdef MODEL400
i_var = (/8,5/) ! Te, rho
#else
i_var = (/6,5/) ! T, rho
#endif
#endif

call interp_PRZ(node_list,element_list,&
  particle%i_elm,i_var,2,particle%x(1),particle%x(2),particle%x(3),&
  P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
#ifdef MODEL400
temperature = max(P(1), 1.d-6) ! fix for negative temperatures
#else
temperature = max(P(1), 1.d-6) * 0.5d0 ! temperature is sum of electron and ion temperatures (assumed equal), unless model400
#endif
!if (P(2) .lt. 0) write(*,*) "WARNING: negative electron density ", P(2), "in update_particle_charge"
P(2) = max(P(2),1.d-3) ! Correct for negative densities

n_e = log10(P(2)*central_density)+20.d0 !log10 of electron number density in m^-3
T_e = log10(temperature/(K_BOLTZ*MU_ZERO*central_density))-20.d0 !log10 of electron temperature in K

particle%q = int(new_charge(int(particle%q,4), ad, n_e, T_e, timestep),1)
end subroutine
end module mod_jorek_ionisation_recombination

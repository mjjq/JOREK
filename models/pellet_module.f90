module pellet_module

use constants
use data_structure
use phys_module


real*8 :: total_pellet_particles   !< the (total) pellet particles added in this timestep
real*8 :: total_plasma_particles   !< the total plasma density (before this timestep)
real*8 :: total_pellet_volume      !< the volume of the simulated pellet in this timestep

real*8 :: phys_pellet_volume       !< the physical pellet radius (in m^3)
real*8 :: pellet_volume            !< approximated value of simulated pellet volume
real*8 :: pellet_atomic            !< atomic number of pellet mass

real*8 :: phys_ablation            !< physical ablation rate (non normalised)



real*8, allocatable  :: xtime_pellet_R(:)
real*8, allocatable  :: xtime_pellet_Z(:)
real*8, allocatable  :: xtime_pellet_psi(:)
real*8, allocatable  :: xtime_pellet_particles(:)
real*8, allocatable  :: xtime_phys_ablation(:)

contains

subroutine pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                          pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                          R,Z,psi,phi, r0, T0, central_density, pellet_particles, pellet_density, pellet_volume, &
                          particle_source, volume_source)

implicit none

!input variables
real*8 :: R, Z, psi, phi            ! position whereh the particle source will be calculated
real*8 :: T0, r0                    ! local temperature and mass density (JOREK normalised)
real*8 :: central_density           !< central plasma density (in units 10^20 m^-3)
real*8 :: pellet_particles          !< total number of particles in the pellet
real*8 :: pellet_density            !< pellet density (units 10^20 m^-3)
real*8 :: pellet_amplitude          !< amplitude of paricle source (when not using ablation model)
real*8 :: pellet_R, pellet_Z        !< position of the pellet (phi=0)
real*8 :: pellet_phi                !< length of pellet in toroidal direction
real*8 :: pellet_radius             !< pellet size (radius) in poloidal plane
real*8 :: pellet_sig, pellet_length !< sigmas of pellet source in poloidal and toroidal direction
real*8 :: pellet_ellipse            !< ellipticity of the pellet source in the poloidal plane
real*8 :: pellet_theta              !< orientation of the pellet ellipse in the poloidal plane
real*8 :: pellet_psi, pellet_delta_psi
real*8 :: pellet_volume

!output variables
real*8 :: particle_source           !< particle source (JOREK normalised units)
real*8 :: volume_source             !< volume of the pellet source (variable used to integrate total pellet volume)

!local variables
real*8  :: radius, atn, atn_psi, atn_phi, atomic_mass, ablation_rate

particle_source = 0.d0
volume_source   = 0.d0

if (pellet_amplitude .gt. 0.) then             ! use the fixed source pellet model 

pellet_particles = 0.0

  if (phi .gt. PI) phi = 2*PI - phi

  radius = sqrt(  (cos(pellet_theta)**2 + 1./pellet_ellipse**2 * sin(pellet_theta)**2)*(R-pellet_R)**2  &
                + (sin(pellet_theta)**2 + 1./pellet_ellipse**2 * cos(pellet_theta)**2)*(Z-pellet_Z)**2  &
                + 2.*(R-pellet_R)*(Z-pellet_Z)*sin(pellet_theta)*cos(pellet_theta) * (1./pellet_ellipse**2 - 1.) )

  atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
  atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
  atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))

  particle_source = pellet_amplitude * atn * atn_phi * atn_psi

!S.F. modified here for introducing moving pellet...
else if (pellet_particles .gt. 0.) then

  if (phi .gt. PI) phi = 2*PI - phi

  radius = sqrt(  (cos(pellet_theta)**2 + 1./pellet_ellipse**2 * sin(pellet_theta)**2)*(R-pellet_R)**2  &
                + (sin(pellet_theta)**2 + 1./pellet_ellipse**2 * cos(pellet_theta)**2)*(Z-pellet_Z)**2  &
                + 2.*(R-pellet_R)*(Z-pellet_Z)*sin(pellet_theta)*cos(pellet_theta) * (1./pellet_ellipse**2 - 1.) )

  atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
  atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))
  atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))

!  pellet_volume = PI * pellet_radius**2 * pellet_R * pellet_phi ! simulated pellet volume

  phys_pellet_volume = pellet_particles /pellet_density         ! physical pellet volume 

! the number of particles ablated from the physical pellet in units 10^20 m^-3 per unit of JOREK time

  pellet_atomic = 2.d0

!----------------- model from Kuteev (see Polevoi, PPCF2008)
 ! ablation_rate = 1.62d5 * central_density**(-0.77) * T0**(1.72) * r0**(0.45) * phys_pellet_volume**(0.48) * pellet_atomic**0.217
!----------------- NGS model Parks (see Gal NF2008)
  ablation_rate = 2.01d4 * central_density**(-0.81) * max(T0,0.d0)**(1.64) * max(r0,0.d0)**(0.33) * phys_pellet_volume**(0.44) * pellet_atomic**0.5

! particle source in JOREK normalisation

  particle_source = ablation_rate / central_density  * atn * atn_phi / pellet_volume

  volume_source   = atn * atn_phi

  if(volume_source .lt. 0.d0) then
!    print*, 'volume_source is negative. volume_source=', volume_source
  endif

  if(ablation_rate .lt. 0.d0) then
!    print*, 'ablation rate is negative. step.'
  endif
  

end if

return
end subroutine pellet_source2

subroutine update_pellet(my_id,node_List,element_list)
!******************************************************************************
! routine updates the pellet position and the size of the simulated           *
! and physical pellet sizes (from the integral of the pellet particle source) *
!******************************************************************************

use constants
use data_structure
use phys_module
use mpi_mod
implicit none


type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: psi_axis, psi_bnd
integer :: my_id, ierr
real*8  :: V_normalisation, density, density_in, density_out, pressure,pressure_in,pressure_out

real*8 :: R_out, Z_out, s_out, t_out, P0_s,P0_t,P0_st,P0_ss,P0_tt
integer :: i_elm, ifail

if (pellet_amplitude .gt. 0) return

call Integrals_3D(my_id, node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO) ! assumes Deuterium!

pellet_R = pellet_R + pellet_velocity_R * tstep / V_normalisation
pellet_Z = pellet_Z + pellet_velocity_Z * tstep / V_normalisation


phys_ablation = total_pellet_particles * central_density / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO)

total_pellet_particles = total_pellet_particles * central_density * tstep 
total_plasma_particles = total_plasma_particles * central_density          ! undo normalisation


call find_RZ(node_list,element_list,pellet_R,pellet_Z,R_out,Z_out,i_elm,s_out,t_out,ifail)
call interp(node_list,element_list,i_elm,1,1,s_out,t_out,pellet_psi,P0_s,P0_t,P0_st,P0_ss,P0_tt)

if (my_id .eq. 0) then

    pellet_particles = max(pellet_particles - total_pellet_particles, 0.d0)

    write(*,'(A,4e14.6)') ' pellet (R,Z) =', pellet_R, pellet_Z,pellet_velocity_R/V_normalisation,pellet_velocity_Z/V_normalisation
    write(*,'(A,6e14.6,A)') ' total particles added in this step : ', pellet_R, pellet_Z,pellet_particles,total_pellet_particles, total_plasma_particles,phys_ablation,' [10^20]'
    write(*,'(A,4e14.6)') ' remaining particles in pellet      : ', pellet_particles
    write(*,'(A,4e14.6)') ' pellet volume (sim,phys)           : ', total_pellet_volume,pellet_particles/pellet_density

else 
  pellet_particles = 0.0
end if
 
call MPI_Bcast(pellet_particles,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

return 
 
end subroutine update_pellet

subroutine update_spi(my_id,node_List,element_list)
!******************************************************************************
! routine updates the shattered pellet position and the size of the simulated *
! and physical pellet sizes (from the integral of the pellet particle source) *
!******************************************************************************

use constants
use data_structure
use phys_module
use mpi_mod
use mgi_module
use corr_neg

implicit none


type (type_node_list)    :: node_list
type (type_element_list) :: element_list

! --- Local variables
integer :: my_id, ierr
real*8  :: V_normalisation, density, density_in, density_out, pressure,pressure_in,pressure_out

real*8  :: R_out, Z_out
integer :: i_elm, ifail

integer :: i

real*8, dimension(3) :: P, P_s, P_t, P_phi
real*8  :: R, R_s, R_t, Z, Z_s, Z_t
real*8  :: s_out,t_out

real*8  :: n_SI, T_eV, n_corr, T_corr, n_imp_SI, ne_SI
real*8  :: t_norm, spi_Vel_totref
real*8  :: spi_delta_phi, spi_Vel_R_tmp, spi_Vel_phi_tmp, spi_phi_inj

!   -Mean impurity ionization state and related quantities
real*8     :: Z_imp, beta_imp, mu_imp

  spi_delta_phi   = 0.
  spi_Vel_R_tmp   = 0.
  spi_Vel_phi_tmp = 0.
  spi_phi_inj     = mgi_phi

  V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO) ! assumes Deuterium!
  t_norm          = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)

  spi_Vel_totref  = sqrt(spi_Vel_Rref**2+spi_Vel_Zref**2+spi_Vel_RxZref**2)
    
  spi_phi_inj     = mgi_phi + mgi_phi_rotate - spi_L_inj * (spi_Vel_RxZref/spi_Vel_totref)/mgi_R

  if (spi_phi_inj >= 2.*PI) then
    spi_phi_inj   = mod(spi_phi_inj,2.*PI)
  else if (spi_phi_inj < 0.) then
    spi_phi_inj   = mod(spi_phi_inj,2.*PI) + 2.*PI
  end if

  do i=1, n_spi

    spi_delta_phi          = pellets(i)%spi_phi - spi_phi_inj
    spi_Vel_R_tmp          = pellets(i)%spi_Vel_R * cos(spi_delta_phi) &
                             + pellets(i)%spi_Vel_RxZ * sin(spi_delta_phi)
    spi_Vel_phi_tmp        = pellets(i)%spi_Vel_RxZ * cos(spi_delta_phi) &
                             - pellets(i)%spi_Vel_R * sin(spi_delta_phi)
    spi_Vel_phi_tmp        = spi_Vel_phi_tmp / pellets(i)%spi_R


    pellets(i)%spi_R       = pellets(i)%spi_R + spi_Vel_R_tmp * tstep / V_normalisation
    pellets(i)%spi_Z       = pellets(i)%spi_Z + pellets(i)%spi_Vel_Z * tstep / V_normalisation
    pellets(i)%spi_phi     = pellets(i)%spi_phi + spi_Vel_phi_tmp * tstep / V_normalisation

    if (toroidal_rotation) then
      pellets(i)%spi_phi     = pellets(i)%spi_phi + tor_frequency * 2. * PI * tstep / V_normalisation
    end if

    pellets(i)%spi_Vel_R   = pellets(i)%spi_Vel_R
    pellets(i)%spi_Vel_Z   = pellets(i)%spi_Vel_Z
    pellets(i)%spi_Vel_RxZ = pellets(i)%spi_Vel_RxZ

    if (pellets(i)%spi_phi >= 2.*PI) then
      pellets(i)%spi_phi   = mod(pellets(i)%spi_phi,2.*PI)
    else if (pellets(i)%spi_phi < 0.) then
      pellets(i)%spi_phi   = mod(pellets(i)%spi_phi,2.*PI) + 2.*PI
    end if

    if (pellets(i)%spi_radius > 0.0) then
      pellets(i)%spi_radius = pellets(i)%spi_radius - t_norm * tstep * &
                              (pellets(i)%spi_abl / (4.d0 * PI * pellets(i)%spi_radius**2.d0 *    &
                              pellet_density * 1.d20))

      if (pellets(i)%spi_radius < 0.d0) then
        pellets(i)%spi_radius = 0.d0
      end if
    end if

    if (my_id == 0.) then
      if (index_now > 1) then
        xtime_spi_ablation(i,index_now) = xtime_spi_ablation(i,index_now-1) + t_norm * tstep * pellets(i)%spi_abl
      else
        xtime_spi_ablation(i,index_now) = t_norm * tstep * pellets(i)%spi_abl
      end if
    end if

    if (flag_spi == 0) then
      pellets(i)%spi_abl   = mgi_amplitude
    elseif (flag_spi >= 1 .and. flag_spi <= 2 ) then

      call find_RZ(node_list,element_list,pellets(i)%spi_R,pellets(i)%spi_Z,&
                   R_out,Z_out,i_elm,s_out,t_out,ifail)

      ! In case the shards are outside of the domain
      if (ifail == 99 .or. ifail == 999) then
        pellets(i)%spi_abl = 0.
        cycle
      else if (ifail /= 0) then
        write(*,*) "Something Wrong in find_RZ!! my_id = ", my_id, i_elm, ifail
        stop
      end if

      call interp_PRZ(node_list,element_list,i_elm,[5,6,8],3,s_out,t_out,pellets(i)%spi_phi,&
                      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)

      ! Now, P(1) represents mass density and P(2) represents temperature, P(3)
      ! is the impurity density
      ! Correct any possible negative values!

      !n_corr         = corr_neg_dens(P(1))
      !T_corr         = corr_neg_temp(P(2))

      ! Reminder, temperature should be divided by 2 since T = T_e + T_i and T_e
      ! = T_i
      !n_SI           = n_corr * 1.d20 * central_density
      !T_eV           = T_corr / (2.d0* EL_CHG * MU_ZERO * central_density * 1.d20)
      
      n_SI           = P(1) * 1.d20 * central_density
      if (n_SI < 0.) n_SI = 0.

      T_eV           = P(2) / (2.d0* EL_CHG * MU_ZERO * central_density * 1.d20)
      if (T_eV < 0.) T_eV = 0.

      ! This is only used for 501 as impurity density, 500 don't use this
      ! variable
      n_imp_SI           = P(3) * 1.d20 * central_density
      if (n_imp_SI < 0.) n_imp_SI = 0.      

      ! NGS model
      if (flag_spi == 1) then
        pellets(i)%spi_abl    = 4.12d16 * (pellets(i)%spi_radius**(4.0/3.0)) * (n_SI**(1.0/3.0)) * &
                               (T_eV**1.64)
        if (my_id == 0 .and. pellets(i)%spi_radius > 0.0 .and. mod(index_now,20)==0) then
          write(*,*) "Check Point, n_SI, T_eV = ", n_SI, T_eV
        end if
      else if (flag_spi == 2) then
        select case ( trim(gas_type) )
          case('D2')
            ne_SI   = n_SI
            ! The scaling law is in gauss unit
            pellets(i)%spi_abl = 3.9d14 * ((pellets(i)%spi_radius*1.d2)**1.455) &
                                 * ((ne_SI*1.d-6)**0.455) * (T_eV**1.679)
          case('Ar')
            if (flag_adas .and. T_eV >= 1.) then
              ! As with element_matrix, mimick density as 1.d20
              call imp_cor(1)%interp(density=20.,temperature=log10(T_eV*EL_CHG/K_BOLTZ),z_out=Z_imp)
            else if (flag_adas .and. T_eV < 1.) then
              Z_imp = 0.
            else
              Z_imp = 10.
            end if
            mu_imp             = 1./20. ! Argon mass = 40 u and main ion (D) mass = 2 u
            beta_imp           = mu_imp*Z_imp - 1.
            ne_SI              = n_SI + beta_imp * n_imp_SI

            if (ne_SI<0.) ne_SI = 0.        
            ! The scaling law is in gauss unit
            pellets(i)%spi_abl = 2.5d13 * ((pellets(i)%spi_radius*1.d2)**1.451) &
                                 * ((ne_SI*1.d-6)**0.451) * (T_eV**1.679)
          ! Using general scaling law of Sergeev for Neon
          case('Ne')
            if (flag_adas .and. T_eV >= 1.) then
              ! As with element_matrix, mimick density as 1.d20
              call imp_cor(1)%interp(density=20.,temperature=log10(T_eV*EL_CHG/K_BOLTZ),z_out=Z_imp)
            else if (flag_adas .and. T_eV < 1.) then
              Z_imp = 0.
            else
              Z_imp = 6.
            end if
            mu_imp             = 1./10. ! Neon mass = 20 u and main ion (D) mass = 2 u
            beta_imp           = mu_imp*Z_imp - 1.
            ne_SI              = n_SI + beta_imp * n_imp_SI

            if (ne_SI<0.) ne_SI = 0.
            ! The scaling law is in gauss unit
            ! The sublimation energy for Ne is 0.02 eV
            pellets(i)%spi_abl = 1.94d14 * ((pellets(i)%spi_radius*1.d2)**1.44) &
                                 * ((ne_SI*1.d-6)**0.45) * (T_eV**1.72)         &
                                 * (0.02**(-0.16)) * (20.**(-0.28))             &
                                 * (10.**(-0.56)) * ((2./3.)**0.28)

          case default
            write(*,*) '!! Gas type "', trim(gas_type), '" unknown !!'
            write(*,*) '=> We assume the gas is D2.'
            pellets(i)%spi_abl = 3.9d14 * ((pellets(i)%spi_radius*1.d2)**1.455) &
                                 * ((n_SI*1.d-6)**0.455) * (T_eV**1.679)
        end select
        if (my_id == 0 .and. pellets(i)%spi_radius > 0.0 .and. mod(index_now,20)==0) then
          write(*,*) "Check Point, ne_SI, T_eV = ", ne_SI, T_eV
        end if
      end if
    else
      pellets(i)%spi_abl    = 0.d0
    end if
   
    if (my_id == 0.) then
      xtime_spi_ablation_rate(i,index_now) = pellets(i)%spi_abl
    end if

  end do

  if (toroidal_rotation) then
    mgi_phi_rotate  = mgi_phi_rotate + tor_frequency * 2. * PI * tstep / V_normalisation
  end if


  !write(*,'(A,4e14.6)') ' pellet (R,Z) =', spi_R, spi_Z,spi_Vel_R/V_normalisation,spi_Vel_Z/V_normalisation
  
  if (my_id == 0 .and. mod(index_now,20) == 0) then

    do i=1, 20 !n_spi
      if (pellets(i)%spi_radius > 0.0) then
        write(*,*) "Pellet number: ", i
        write(*,*) "Pellet coordinates (R,Z,phi) = ", pellets(i)%spi_R, pellets(i)%spi_Z, pellets(i)%spi_phi
        write(*,*) "Pellet velocity (R,Z,phi) = ", pellets(i)%spi_Vel_R, pellets(i)%spi_Vel_Z, &
                                                   pellets(i)%spi_Vel_RxZ
        write(*,*) "Pellet ablation (radius,abl) = ", pellets(i)%spi_radius, pellets(i)%spi_abl
      end if
    end do
  end if

  !call Integrals_3D(my_id, node_list,element_list,density,&
        !density_in,density_out,pressure,pressure_in,pressure_out)

return 
 
end subroutine update_spi

!> This function creates a derived MPI type for the pellets and returns it (in honor of Daan)
!! If it already exists the old handle is returned
function get_pellet_derived_type() result(dtype_out)
  use mpi_mod
  use mod_parameters

  implicit none

  integer               :: ierr, dtype_out
  integer, save         :: dtype
  logical, save         :: dtype_set = .false.

  integer :: len(8) = (/1,1,1,1,1,1,1,1/), t(8) = (/ &
    MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8, &
    MPI_REAL8,MPI_REAL8,MPI_REAL8/) ! MPI_INTEGER1 == MPI_LOGICAL1

  integer(kind=MPI_ADDRESS_KIND) :: base, disp(8)
  type(type_SPI) :: sample_pellet

  dtype_out = dtype
  if (dtype_set) return

  ! Get memory addresses in the type
  call MPI_Get_address(sample_pellet,             base,    ierr)
  call MPI_Get_address(sample_pellet%spi_R,       disp(1), ierr)
  call MPI_Get_address(sample_pellet%spi_Z,       disp(2), ierr)
  call MPI_Get_address(sample_pellet%spi_phi,     disp(3), ierr)
  call MPI_Get_address(sample_pellet%spi_Vel_R,   disp(4), ierr)
  call MPI_Get_address(sample_pellet%spi_Vel_Z,   disp(5), ierr)
  call MPI_Get_address(sample_pellet%spi_Vel_RxZ, disp(6), ierr)
  call MPI_Get_address(sample_pellet%spi_radius,  disp(7), ierr)
  call MPI_Get_address(sample_pellet%spi_abl,     disp(8), ierr)

  ! Rebase to particle memory beginning
  disp = disp - base

  ! Commit the structured type
  call MPI_Type_create_struct(8, len, disp, t, dtype, ierr)
  call MPI_Type_commit(dtype, ierr)

  ! Set the save bit
  dtype_set = .true.
  dtype_out = dtype
  return
end function get_pellet_derived_type

subroutine pellet_source(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                         pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
			 R,Z,psi,phi,particle_source)

use constants

implicit none

real*8 :: R, Z, psi, phi, particle_source
real*8 :: pellet_amplitude, pellet_R, pellet_Z, pellet_phi, pellet_radius, pellet_sig, pellet_length
real*8 :: pellet_psi, pellet_delta_psi, radius, atn, atn_psi, atn_phi

if (phi .gt. PI) phi = 2*PI - phi

radius = sqrt((R-pellet_R)**2 + (Z-pellet_Z)**2)

atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))

atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))

atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))

particle_source = pellet_amplitude * atn * atn_phi * atn_psi

return
end subroutine pellet_source
end module pellet_module

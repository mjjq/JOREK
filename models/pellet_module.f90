module pellet_module

  use mod_interp
  use constants
  use data_structure

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
#if _OPENMP >= 201511
  !$omp declare simd uniform(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
  !$omp                           pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
  !$omp                          R,Z, r0, T0, central_density, pellet_particles, pellet_density, pellet_volume)
#endif
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

  !> Update the pellet position, size and ablation rate
  subroutine update_pellet(my_id,node_List,element_list)
  
    use constants
    use data_structure
    use phys_module
    use mpi_mod
    implicit none  
    
    type (type_node_list), intent(in)    :: node_list
    type (type_element_list), intent(in) :: element_list
    
    real*8              :: psi_axis, psi_bnd
    integer, intent(in) :: my_id
    integer             :: ierr, i_elm, ifail
    real*8              :: V_normalisation, density, density_in, density_out, pressure, pressure_in, pressure_out   
    real*8              :: R_out, Z_out, s_out, t_out, P0_s, P0_t, P0_st, P0_ss, P0_tt
    
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

  !> Update the shattered pellet position, size and ablation rate
  subroutine update_spi(my_id,node_List,element_list)
  
    use constants
    use data_structure
    use phys_module
    use mpi_mod
#if (JOREK_MODEL == 500 || JOREK_MODEL == 501 || JOREK_MODEL == 555)
      use mod_neutral_source
#endif
    use corr_neg
    
    implicit none
    
    
    type (type_node_list), intent(in)    :: node_list
    type (type_element_list), intent(in) :: element_list
    
    ! --- Local variables
    real*8               :: psi_axis, psi_bnd
    integer, intent(in)  :: my_id
    integer              :: ierr, i_elm, ifail, i
    real*8               :: V_normalisation, density, density_in, density_out, pressure, pressure_in, pressure_out    
    real*8               :: R_out, Z_out, R, R_s, R_t, Z, Z_s, Z_t, s_out, t_out    
    real*8, dimension(2) :: P, P_s, P_t, P_phi   
    real*8               :: n_SI, T_eV, n_corr, T_corr, t_norm, spi_Vel_totref
    ! Temp variables for SPI, used to advance the R, Z, phi position of pellets for
    ! given R, Z, RxZ velocity and a calculated apex of the trajectory spreading cone 
    real*8               :: spi_delta_phi, spi_Vel_R_tmp, spi_Vel_phi_tmp, spi_phi_inj
  
    spi_delta_phi   = 0.
    spi_Vel_R_tmp   = 0.
    spi_Vel_phi_tmp = 0.
    spi_phi_inj     = ns_phi
    
    V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO)
    t_norm          = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
    
    spi_Vel_totref  = sqrt(spi_Vel_Rref**2+spi_Vel_Zref**2+spi_Vel_RxZref**2)
        
    spi_phi_inj     = ns_phi + ns_phi_rotate - spi_L_inj * (spi_Vel_RxZref/spi_Vel_totref)/ns_R
    
    if (spi_phi_inj >= 2.*PI) then
      spi_phi_inj   = mod(spi_phi_inj,2.*PI)
    else if (spi_phi_inj < 0.) then
      spi_phi_inj   = mod(spi_phi_inj,2.*PI) + 2.*PI
    end if
    
    loop_over_shards:do i = 1, n_spi

      ! Update position     
      spi_delta_phi          = pellets(i)%spi_phi - spi_phi_inj
      spi_Vel_R_tmp          = pellets(i)%spi_Vel_R * cos(spi_delta_phi) &
                               + pellets(i)%spi_Vel_RxZ * sin(spi_delta_phi)
      spi_Vel_phi_tmp        = pellets(i)%spi_Vel_RxZ * cos(spi_delta_phi) &
                               - pellets(i)%spi_Vel_R * sin(spi_delta_phi)
      spi_Vel_phi_tmp        = spi_Vel_phi_tmp / pellets(i)%spi_R    
    
      pellets(i)%spi_R       = pellets(i)%spi_R + spi_Vel_R_tmp * tstep / V_normalisation
      pellets(i)%spi_Z       = pellets(i)%spi_Z + pellets(i)%spi_Vel_Z * tstep / V_normalisation
      pellets(i)%spi_phi     = pellets(i)%spi_phi + spi_Vel_phi_tmp * tstep / V_normalisation
    
      if (spi_tor_rot) then
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

      ! Update size    
      if (my_id == 0) then
        if (index_now > 1) then
          xtime_spi_ablation(i,index_now) = xtime_spi_ablation(i,index_now-1) + t_norm * tstep * pellets(i)%spi_abl
        else
          xtime_spi_ablation(i,index_now) = t_norm * tstep * pellets(i)%spi_abl
        end if
      end if

      ! Update ablation rate    
      if (spi_abl_model == 0) then ! Constant ablation rate
        pellets(i)%spi_abl = ns_amplitude
      elseif (spi_abl_model >= 1 .and. spi_abl_model <= 2) then ! NGS models (see Wiki: https://www.jorek.eu/wiki/doku.php?id=spi_tutorial)

        ! Get local n_e, T_e    
        call find_RZ(node_list,element_list,pellets(i)%spi_R,pellets(i)%spi_Z,&
                     R_out,Z_out,i_elm,s_out,t_out,ifail)
    
        if (ifail == 99 .or. ifail == 999) then ! The shard is outside of the domain
          pellets(i)%spi_abl = 0.
          cycle
        else if (ifail /= 0) then
          write(*,*) "Something wrong in find_RZ!! my_id = ", my_id, i_elm, ifail
          stop
        end if
    
        call interp_PRZ(node_list,element_list,i_elm,[5,6],2,s_out,t_out,pellets(i)%spi_phi,&
                        P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
          
        n_SI = P(1) * 1.d20 * central_density
        n_SI = max(n_SI,0.d0)

        ! Reminder, temperature should be divided by 2 since T = T_e + T_i and T_e = T_i    
        T_eV = P(2) / (2.d0* EL_CHG * MU_ZERO * central_density * 1.d20) 
        T_eV = max(T_eV,0.d0)
        
        if (my_id == 0 .and. pellets(i)%spi_radius > 0.0 .and. mod(index_now,20)==0) then
          write(*,*) "Check Point, n_SI, T_eV = ", n_SI, T_eV
        end if
          
        if (spi_abl_model == 1) then ! NGS model version 1 (K. Gal et al., see http://iopscience.iop.org/article/10.1088/0029-5515/48/8/085005?pageTitle=IOPscience)
          pellets(i)%spi_abl = 4.12d16 * (pellets(i)%spi_radius**(4.0/3.0)) * (n_SI**(1.0/3.0)) * (T_eV**1.64)
        else if (spi_abl_model == 2) then ! NGS model version 2 (V. Sergeev et al., see https://link.springer.com/article/10.1134/S1063780X06050023)
          pellets(i)%spi_abl = 3.9d14 * ((pellets(i)%spi_radius*1.d2)**(1.455)) * ((n_SI*1.d-6)**(0.455)) * (T_eV**1.679)
        end if
      else
        write(*,*) "Forbidden value for spi_abl_model: ", spi_abl_model
        write(*,*) "=> EXIT"
        stop
      end if
       
      if (my_id == 0) then
        xtime_spi_ablation_rate(i,index_now) = pellets(i)%spi_abl
      end if
    
    end do loop_over_shards
    
    if (spi_tor_rot) then
      ns_phi_rotate  = ns_phi_rotate + tor_frequency * 2. * PI * tstep / V_normalisation
    end if
    
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
    
    return 
   
  end subroutine update_spi

  !> Initializes the shattered pellet position, velocity and size
  subroutine init_spi(my_id)
  
    use constants
    use data_structure
    use phys_module
    use mpi_mod
#if (JOREK_MODEL == 500 || JOREK_MODEL == 501 || JOREK_MODEL == 555)
    use mod_neutral_source
#endif
    use corr_neg
    
    implicit none
    
    integer, intent(in) :: my_id
    integer             :: ierr,err,i
    
    logical             :: ferr
    
    real*8  :: n_SI, T_eV, n_corr, T_corr
    real*8  :: spi_gd_angle_01, spi_gd_angle_02        !The dispersion angles for each spi
    real*8  :: spi_rotation_01, spi_rotation_02        !The rotation angle from spi coordinate to real coordinate
    real*8  :: spi_Vel_totref, spi_Vel_i
    real*8  :: spi_Vel_x, spi_Vel_y, spi_Vel_z         !Spi velocity in injection coordinate
    real*8  :: spi_R_inj, spi_Z_inj, spi_phi_inj       !Representing the shattering point of the pellet
                                                       !The apex of the spreading cone
    real*8  :: sign_corr
    real*8, allocatable :: rnd(:)                      !The random number array 
    real*8, allocatable :: shard_size(:)               !The shard size array

    real*8  :: size_beta                               ! The characteristic shard size    
    real*8  :: V_shard_norm                            ! The normalized (by size_beta) volume of shards
    
    if (allocated(pellets)) then
      deallocate(pellets)
    end if
  
    if (n_spi >= 1) then
 
      allocate (pellets(n_spi))  ! Dynamically allocate memeries for pellets
 
      ! Read normalized shard size distribution (if given in file) and calculate shard radius normalization factor size_beta
      if (spi_shard_file /= 'none') then 
  
        if (allocated(shard_size)) then
          deallocate(shard_size)
        end if
        allocate (shard_size(n_spi))  ! Dynamically allocate memeries for shard sizes
        shard_size(:) = 0.0
  
        size_beta    = 0.0
        V_shard_norm = 0.0
        inquire(file=trim(spi_shard_file), exist=ferr) ! Check if the file exists
        if (ferr) then
          open(42,file=trim(spi_shard_file),status="OLD",action="READ")
          read(42,*)  shard_size(1:n_spi)
          close(42)
        else
          write(*,*) "WARNING!!! Shard size file does not exist!"
          if (index_now == 0) then
            deallocate(shard_size)
            deallocate(pellets)
            stop
          else
            write(*,*) "The shard information will be overwritten by restart read anyway, proceed assuming no shard size file is given."
            spi_shard_file = 'none'
            shard_size(:) = 1.0
          end if
        end if

        do i = 1, n_spi
          V_shard_norm = V_shard_norm + (4./3.) * PI * (shard_size(i)**3)
        end do

        size_beta    = (spi_quantity / (V_shard_norm)*pellet_density*1.d20) ** (-1./3.)
        write(*,*) "Characteristic shard size (m):", size_beta

      end if

      ! Initialize shard radius
      do i = 1, n_spi 
        if (spi_shard_file == 'none') then
          pellets(i)%spi_radius = (spi_quantity / (n_spi*(4.*PI/3.)*pellet_density*1.d20))**(1./3.)
        else
          pellets(i)%spi_radius = shard_size(i)/size_beta
        end if
      end do

      ! Initialize shard velocity and position    
  
      ! Determine rotation angles from SPI coordinate system x, y, z to R, Z, RxZ
      ! with the reference direction of spi injection being the z axis, while y axis locates within the 
      ! same surface as Z and z. The rotational transform from x, y, z to R, Z, RxZ is
      ! as the following: first, we rotate the system around x axis clockwise, facing
      ! the positive x direction, for spi_rotation_01 to get coordinate X', Y', Z'. 
      ! Then we further rotate around Y' clockwise, facing the positive Y' direction for
      ! spi_rotation_02 to acquire R, Z, RxZ. Hence we have:
      ! R   = cos(spi_rotation_02)*x - sin(spi_rotation_02)*(-sin(spi_rotation_01)*y + cos(spi_rotation_01)*z)
      ! Z   = cos(spi_rotation_01)*y + sin(spi_rotation_01)*z
      ! RxZ = sin(spi_rotation_02)*x + cos(spi_rotation_02)*(-sin(spi_rotation_01)*y + cos(spi_rotation_01)*z)  
      spi_Vel_totref  = sqrt(spi_Vel_Rref**2+spi_Vel_Zref**2+spi_Vel_RxZref**2)  
      spi_R_inj       = ns_R - spi_L_inj * (spi_Vel_Rref/spi_Vel_totref)
      spi_Z_inj       = ns_Z - spi_L_inj * (spi_Vel_Zref/spi_Vel_totref)
      spi_phi_inj     = ns_phi - spi_L_inj * (spi_Vel_RxZref/spi_Vel_totref)/ns_R  
      spi_rotation_01 = asin(spi_Vel_Zref/spi_Vel_totref)
      if (cos(spi_rotation_01) == 0.) then
        spi_rotation_02 = 0.
      else
        spi_rotation_02 = acos(spi_Vel_RxZref/(spi_Vel_totref*cos(spi_rotation_01)))
      end if 
      write(*,*) "Rotation angles from SPI coordinate system to R, Z, RxZ: ", spi_rotation_01, spi_rotation_02

      ! Generate a random number array rnd that contains two random angles representing the
      ! velocity direction spread, and one the random speed. Those random numbers uniquely
      ! define a random velocity of the shard, which is then transformed into the
      ! R, Z, RxZ space.
      if (allocated(rnd)) deallocate(rnd)
      allocate (rnd(3*n_spi))  ! Dynamically allocate memeries for randoms   
      CALL random_number(rnd)  
    
      spi_Vel_diff = abs(spi_Vel_diff) ! To avoid negative velocity spread
  
      do i = 1, n_spi
  
        spi_gd_angle_01 = rnd(3 * i - 2) * spi_angle / 2.0
        spi_gd_angle_02 = rnd(3 * i - 1) * 2. * PI
        spi_Vel_i       = (rnd(3*i)-0.5) * spi_Vel_diff + spi_Vel_totref 
  
        !write(*,*) "Random angle:", i, spi_gd_angle_01, spi_gd_angle_02
        spi_Vel_x       = spi_Vel_i * sin(spi_gd_angle_01) * cos(spi_gd_angle_02)
        spi_Vel_y       = spi_Vel_i * sin(spi_gd_angle_01) * sin(spi_gd_angle_02)
        spi_Vel_z       = spi_Vel_i * cos(spi_gd_angle_01)
  
        pellets(i)%spi_Vel_R   = spi_Vel_x * cos(spi_rotation_02) &
                                 - sin(spi_rotation_02) * (-sin(spi_rotation_01)*spi_Vel_y &
                                 + cos(spi_rotation_01)*spi_Vel_z)
        pellets(i)%spi_Vel_Z   = cos(spi_rotation_01) * spi_Vel_y &
                                 + sin(spi_rotation_01) * spi_Vel_z
        pellets(i)%spi_Vel_RxZ = spi_Vel_x * sin(spi_rotation_02) &
                                 - cos(spi_rotation_02) * (-sin(spi_rotation_01)*spi_Vel_y &
                                 + cos(spi_rotation_01)*spi_Vel_z)
  
        pellets(i)%spi_R       = spi_R_inj + spi_L_inj * (pellets(i)%spi_Vel_R/spi_Vel_totref)
        pellets(i)%spi_Z       = spi_Z_inj + spi_L_inj * (pellets(i)%spi_Vel_Z/spi_Vel_totref)
        pellets(i)%spi_phi     = spi_phi_inj + spi_L_inj * (pellets(i)%spi_Vel_RxZ/spi_Vel_totref)/ns_R
  
        pellets(i)%spi_abl     = 0.0
  
        write(*,'(A,I5,5ES10.2)') ' *** SHATTERED PELLET PARAMETERS :',i, pellets(i)%spi_R, pellets(i)%spi_Z, &
                                pellets(i)%spi_Vel_R, pellets(i)%spi_Vel_Z, pellets(i)%spi_radius
        
      end do

      if (allocated(rnd)) deallocate(rnd)
      if (allocated(shard_size)) deallocate(shard_size)
  
      if (allocated(xtime_spi_ablation)) call tr_deallocate(xtime_spi_ablation,"xtime_spi_ablation",CAT_GRID)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation,1,n_spi,1,nstep,"xtime_spi_ablation")
  
      if (allocated(xtime_spi_ablation_rate)) &
      call tr_deallocate(xtime_spi_ablation_rate,"xtime_spi_ablation_rate",CAT_GRID)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_rate,1,n_spi,1,nstep,"xtime_spi_ablation_rate")
      
    else ! In the illigal case of n_spi < 1
      write(*,*) "......Seriously!? Using the spi flag while n_spi is set to be smaller than 1, double check the input file."
      stop
    end if
    
    return
  end subroutine init_spi


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

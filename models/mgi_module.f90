module mgi_module

  use constants

  real*8 :: total_n_particles_inj
  real*8 :: total_n_particles
  real*8 :: total_n_particles_inj_all

  contains 


  integer function factorial(n)

    implicit none

    integer, intent(in) :: n 
    integer             :: i, Ans

    Ans = 1
 
    do i=1,n
      Ans = Ans * i
    enddo

    factorial = Ans

  end function factorial



  subroutine mgi_source(mgi_amplitude,mgi_R,mgi_Z,mgi_phi,mgi_radius,mgi_sig,mgi_deltaphi,mgi_tor_norm,  &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_mgi,L_tube,R,Z,phi,rhon_source,t_now,                  &
                        JET_MGI,ASDEX_MGI,central_density,central_mass)

  !=================================================================================
  !  This subroutine computes the neutral density source for a realistic Deuterium
  !  MGI in JET (if mgi_timedependent is .t.).
  !  If mgi_timedependent is .f., this routine computes a constant source in time
  !  where the main parameter is mgi_amplitude
  !  More details in the JOREK wiki or by asking A.Fil or E.Nardon
  !=================================================================================

    use phys_module, only: gas_type

    implicit none

    real*8 :: c0_gas                   ! Sound velocity of gas in reservoir
    integer:: n_gas                    ! = 2/(gamma-1) where gamma = heat capacity ratio of gas
    real*8 :: A_gas                    ! Atomic number of gas particles
    real*8 :: mass_gas                 ! Mass of a gas particles
    real*8 :: mol_atom                 ! Number of atoms in a molecular
    real*8 :: radius
    real*8 :: mgi_tor_shape
    real*8 :: mgi_pol_shape
    real*8 :: dphi
    real*8 :: V_mgi
    real*8 :: f_Nbar
    real*8 :: f_dNbar_dt
    real*8 :: mgi_dNinj_dt
    real*8 :: mgi_drhon_dt
    real*8 :: t_loc
    real*8 :: t_norm
    real*8 :: prof_temp
    real*8 :: R_Asdex
    real*8 :: mnum
    real*8 :: kst
    real*8 :: yy
    real*8 :: gam
    real*8 :: dt_open
    real*8 :: N_barlitre
    integer:: k
    real*8, intent(in)  :: R
    real*8, intent(in)  :: Z
    real*8, intent(in)  :: phi
    real*8, intent(in)  :: A_Dmv
    real*8, intent(in)  :: K_Dmv
    real*8, intent(in)  :: V_Dmv
    real*8, intent(in)  :: P_Dmv
    real*8, intent(in)  :: t_now
    real*8, intent(in)  :: t_mgi
    real*8, intent(in)  :: mgi_amplitude
    real*8, intent(in)  :: mgi_R
    real*8, intent(in)  :: mgi_Z
    real*8, intent(in)  :: mgi_phi
    real*8, intent(in)  :: mgi_radius
    real*8, intent(in)  :: mgi_sig
    real*8, intent(in)  :: mgi_deltaphi
    real*8, intent(in)  :: L_tube
    real*8, intent(in)  :: central_density
    real*8, intent(in)  :: central_mass
    real*8              :: DMV_inj_frac
    logical, intent(in) :: JET_MGI
    logical, intent(in) :: ASDEX_MGI
    real*8, intent(out) :: rhon_source
    real*8, intent(in)  :: mgi_tor_norm

    select case ( trim(gas_type) )
      case('D2')
        n_gas  = 5
        A_gas  = 4.
        mol_atom = 2.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(7.d0/5.d0))
      case('Ar')
        n_gas  = 3
        A_gas  = 40.
        mol_atom = 1.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(5.d0/3.d0))
      case('Ne')
        n_gas  = 3
        A_gas  = 20.
        mol_atom = 1.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(5.d0/3.d0))
      case default
        write(*,*) '!! Gas type "', trim(gas_type), '" unknown (in mgi_source.f90) !!'
        write(*,*) '=> We assume the gas is D2.'
        n_gas  = 5
        A_gas  = 4.
        mol_atom = 2.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(7.d0/5.d0))
    end select

    ! ===================================================================
    ! Parameters related to the spatial distribution of the gas source:

    ! A gaussian shape is chosen poloidally
    radius = sqrt((R-mgi_R)**2 + (Z-mgi_Z)**2)
    mgi_pol_shape = exp(-(radius/mgi_radius)**2.d0)  

    ! A gaussian shape is chosen toroidally
    dphi = abs(phi - mgi_phi)
    if (dphi .gt. PI) dphi = 2*PI - dphi  
    mgi_tor_shape = exp(-(dphi/mgi_deltaphi)**2.d0)

    ! Volume used for normalization, which corresponds to the integration in space 
    ! of the product of the above shape functions
    V_mgi  = PI * mgi_R * mgi_tor_norm * mgi_radius**2.d0
    ! ===================================================================

   !==================================================================================================
   ! A shifted time is used in order to start injected gas as soon as t_now = t_mgi 
   ! (note: L_tube/3c0 is the time needed for the gas to propagate in the injection tube).
    t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
    t_loc = (t_now-t_mgi) * t_norm + L_tube/(3.d0 * c0_gas)
   !==================================================================================================

    if (t_loc .gt. 0.) then

      if (JET_MGI) then

       !==================================================================================================
       ! We use here the formulae derived from the eq(8) in the paper of S.A. Bozhenkov - NF 51 (2011) 
       ! which gives the normalized number of particles injected at the exit of the DMV injection tube
       ! as a function of time.
       ! The parameters used are realistic:
       ! A_Dmv: cross sectional area of the injection pipe
       ! K_Dmv: Experimental correction factor to account for the gas expansion close to the tube orifice
       ! L_tube: DMV vacuum injection tube length
       ! V_Dmv: Volume of the DMV reservoir
       ! P_Dmv: Initial pressure in the DMV reservoir, directly linked to the total number of particles
       ! in the reservoir. Expressed in bar here as it is in all MGI experiments. 
       !==================================================================================================

        f_Nbar = 0.d0
        f_dNbar_dt = 0.d0

       ! Calculation of the normalized number of particles injected per unit time, following Bozhenkov
        do k = 0,n_gas+1
          f_Nbar     = f_Nbar + (-1.d0)**(k-1)*factorial(n_gas+1)/(factorial(n_gas+1-k)*factorial(k))*(1-(n_gas*c0_gas*t_loc/L_tube)**(1-k))

          f_dNbar_dt = f_dNbar_dt &
                        + (-1.d0)**(k-1)*factorial(n_gas+1)/(factorial(n_gas+1-k)*factorial(k))*(k-1)*(n_gas*c0_gas*(L_tube)**(-1.d0))**(1-k) &
                          *t_loc**(-k)
        end do

        f_Nbar     = ((1.*n_gas)**n_gas) * ((1.*(n_gas+1.))**(-n_gas-1)) * f_Nbar
        f_dNbar_dt = ((1.*n_gas)**n_gas) * ((1.*(n_gas+1.))**(-n_gas-1)) * f_dNbar_dt

        DMV_inj_frac = A_Dmv * L_tube * K_Dmv * f_Nbar/(V_Dmv)

       ! The gas injection is stopped when the initial number of particles in the reservoir is reached 
       ! if (DMV_inj_frac .gt. 1.d0) then      
       !   f_dNbar_dt = 0.d0
       ! endif

        ! Number of injected particles per unit time, normalized to reservoir content:
        mgi_dNinj_dt = A_Dmv * K_Dmv * L_tube / V_Dmv * f_dNbar_dt

        ! Mass density injected per unit time (SI units):
        mgi_drhon_dt = mgi_dNinj_dt * (P_Dmv * 1.d5/(K_BOLTZ * 293)) * V_Dmv * mass_gas
    
        ! Distribute gas source in space
        rhon_source = mgi_drhon_dt * mgi_pol_shape * mgi_tor_shape / V_mgi

        ! Apply JOREK normalization
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * rhon_source

      elseif (ASDEX_MGI) then

        N_barlitre = (6.02d23*1.d5*1.d-3)/(8.3144d0*293d0)

        R_Asdex = 8314.4d0
        mnum = 20.2d0
        kst = 1.666d0 
        ! kst= 5/3 for noble gas as Ne;
        ! A_Dmv = PI*0.7*0.7*1d-4
        ! V_Dmv = 80.0d-6

        yy = (2.d0/(1+kst))**(1.d0/(kst-1))*(2.d0*kst/(kst+1)*R_Asdex*293.d0/mnum)**0.5d0
    
        gam = A_Dmv/V_Dmv*yy
    
        dt_open = 1.5d-3

        if (t_loc .lt. dt_open) then

          prof_temp = - exp(-t_loc*t_loc/2.d0/dt_open*gam)

          mgi_dNinj_dt = - prof_temp*t_loc*V_Dmv*1.d3*gam*P_Dmv*N_barlitre/dt_open ! Number of injected particles per unit time (not normalised)

        else

          prof_temp = - exp(-(t_loc-dt_open)*gam)*exp(-dt_open/(2*gam))

          mgi_dNinj_dt = - prof_temp*gam*V_Dmv*1.d3*P_Dmv*N_barlitre ! Number of injected particles per unit time (not normalised)
    
        endif

        mgi_drhon_dt =  mgi_dNinj_dt * mass_gas ! Mass density injected per unit time
    
        ! Inverse of the number of particles still in the reservoir, formulae given by G. Pautasso (ASDEX-U)

        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0)*mgi_drhon_dt * mgi_pol_shape  * mgi_tor_shape / V_mgi

      else 

        rhon_source = rhon_source + (mgi_amplitude * mgi_pol_shape * mgi_tor_shape * t_norm &
                      * mass_gas /  (V_mgi * mol_atom * 1.d20 * central_density * central_mass * MASS_PROTON))  
 
      endif

    else

      rhon_source = 0.

    endif

  if (rhon_source < 0.) then
    rhon_source = 0.
  end if


  return
  end subroutine mgi_source


  subroutine update_mgi(my_id,node_list,element_list)

  !=================================================================================================
  ! This routine is used to calculate the total number of neutral particles injected (in nb part/s)
  ! from the start of the simulation and for each timestep.
  ! It also calculate to total number of neutral particles in the plasma.
  !=================================================================================================
    use data_structure
    use phys_module
    use mpi_mod

    implicit none

    type (type_node_list), intent(in)    :: node_list
    type (type_element_list), intent(in) :: element_list

    integer, intent(in) :: my_id 
    integer :: ierr
    real*8  :: density, density_in, density_out, pressure, pressure_in,pressure_out

    call Integrals_3D(my_id, node_list, element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

    total_n_particles_inj_all = total_n_particles_inj_all + total_n_particles_inj*tstep*sqrt(MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)

    if (my_id .eq. 0) then
      write(*,*) 'total neutrals particles injected per second =', total_n_particles_inj
      write(*,*) 'total neutrals particles into the plasma =', total_n_particles
      write(*,*) 'total neutrals particles injected since the start of the simulation = ', total_n_particles_inj_all
      !write(*,*) 'Check of density conservation'
    endif

  end subroutine update_mgi

  subroutine init_imp_adas(my_id)

    use phys_module
    use mod_openadas
    use mod_coronal

    implicit none

    integer, intent(in) :: my_id
    integer             :: err_alloc, i

    character(len=512)  :: adas_suffix     !The suffix of adas data file to be read

    ! Temporary variable for charge state distribution
    integer             :: i_T, i_ion
    real*8, allocatable :: dP_imp_dT(:), P_imp(:)
    real*8              :: T_rad(101) = 0
    real*8              :: Z_imp

    n_adas = 1 ! For now we only trace one species, in the future probably more 

    if (allocated(imp_adas)) then
      deallocate(imp_adas)
    end if

    allocate (imp_adas(n_adas),stat=err_alloc)  !< Dynamically allocate memeries for adas data

    if (err_alloc /= 0) then
      write(*,*) "Error when trying to dynamically allocate memeries for adas data.", my_id
      stop
    else
      if (allocated(imp_cor)) then
        deallocate(imp_cor)
      end if

      allocate (imp_cor(n_adas),stat=err_alloc)  !< Dynamically allocate memeries for adas data
      if (err_alloc /= 0) then
        write(*,*) "Error when trying to dynamically allocate memeries for CE vector.", my_id
        deallocate(imp_adas)
        stop
      else
        do i=1, n_adas
          select case ( trim(gas_type) )
            case('D2')
              write(*,*) "Deuterium adas calculation unsupported for now, disable flag_adas."
              adas_suffix = 'none'
              deallocate(imp_cor)
              deallocate(imp_adas)
              stop
            case('Ar')
              adas_suffix = '89_ar'
            case('Ne')
              adas_suffix = '96_ne'
            case default
              write(*,*) "Unrecognized species, disable flag_adas."
              adas_suffix = 'none'
              deallocate(imp_cor)
              deallocate(imp_adas)
              stop
          end select

          imp_adas(i) = read_adf11(trim(adas_suffix),trim(adas_dir))
          imp_cor(i)  = coronal(imp_adas(i))

          
          ! This is to output a coronal equilibrium charge distribution as a
          ! function of temperature assuming constant density
          open(20,file="charge_distribution.dat")

          write(20,'(A11)') 'temperature (log10(K))', 'charge states'


          do i_T = 0, 100

            T_rad(i_T) = 2. + 0.06*real(i_T,8)

            if (allocated(P_imp)) deallocate(P_imp)
            if (allocated(dP_imp_dT)) deallocate(dP_imp_dT)

            allocate(P_imp(0:imp_adas(i)%n_Z))
            allocate(dP_imp_dT(0:imp_adas(i)%n_Z))
            call imp_cor(1)%interp(density=20.,temperature=T_rad(i_T),p_out=P_imp,z_eff=Z_imp)
            write(20,'(f12.3)',advance='no'), T_rad(i_T)
            do i_ion = 0, imp_adas(i)%n_Z
              write(20,'(f12.5)',advance='no'), P_imp(i_ion)
            end do
            write(20,'(f12.5)'), sum(P_imp)
          end do

          close(20)
        end do
      end if

    end if

  end subroutine init_imp_adas

  subroutine radiation_function(ad, cor, density, temperature, Lrad, dLrad_dTe, dLrad_dNe)

    use phys_module
    use mod_openadas
    use mod_coronal

    implicit none

    type(coronal), intent(in)   :: cor
    type(adf11_all), intent(in) :: ad
    real*8, intent(in)          :: density !< log10 density in m^-3
    real*8, intent(in)          :: temperature !< log10 electron temperature in K

    real*8, intent(out)         :: Lrad ! value of radiation functioni
    real*8, intent(out), optional :: dLrad_dTe, dLrad_dNe ! derivatives of radiation functioni

    real*8, dimension(ad%n_Z)   :: rad     ! The value of radiation function for each charge state
    real*8, dimension(ad%n_Z)   :: dlograd_dlogTe, dlograd_dlogNe ! The loglog gradient of radiation function 
                                                                  ! for each charge state in K and m^-3
    real*8, dimension(ad%n_Z)   :: drad_dTe, drad_dNe             ! The gradient of radiation function
    real*8, dimension(0:cor%n_Z):: p          !< charge state distribution
    real*8, dimension(0:cor%n_Z):: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
    integer :: iz

    dlograd_dlogTe       = ad%PRB%interp_grad_T(density, temperature) + &
                           ad%PLT%interp_grad_T(density, temperature)
    dlograd_dlogNe       = ad%PRB%interp_grad_n(density, temperature) + &
                           ad%PLT%interp_grad_n(density, temperature)

    do iz=1,ad%n_Z
      rad(iz)            = ad%PRB%interp(iz, density, temperature) + &
                           ad%PLT%interp(iz, density, temperature)
      drad_dTe(iz)       = dlograd_dlogTe(iz) * rad(iz) / (10.0**temperature) !log10 terms cancel
      drad_dNe(iz)       = dlograd_dlogNe(iz) * rad(iz) / (10.0**density)     !log10 terms cancel
    enddo ! radiation emitted by atoms at level iz


    Lrad = L2Dinterp(cor%density,cor%temperature,cor%Prad(:,:),density,temperature)
    Lrad = Lrad / (10.0**density) ! This is to recover the radiation coefficient
    
    if (present(dLrad_dTe) .or. present(dLrad_dNe)) then
      call cor%interp(density,temperature,p)
      call cor%interp_gradients(density,temperature,p_Te,p_Ne) 
      if (present(dLrad_dTe)) then
        dLrad_dTe = dot_product(p_Te(1:ad%n_Z),rad) + dot_product(p(1:ad%n_Z),drad_dTe)
      endif
      if (present(dLrad_dNe)) then
        dLrad_dNe = dot_product(p_Ne(1:ad%n_Z),rad) + dot_product(p(1:ad%n_Z),drad_dNe)
      endif


    endif


  end subroutine radiation_function

end module mgi_module

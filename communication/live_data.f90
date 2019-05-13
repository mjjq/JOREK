!> The module contains routines which write certain data toa text file while the code is running.
!!
!! - The file <b>macroscopic_vars.dat</b> is created during the code run and filled with
!!   information about certain run parameters, energy timetraces, growth rates, etc.
!! - The input parameter phys_module::produce_live_data allows to switch the functionality of
!!   this module on (default) or off.
!! - The script extract_live_data.sh in the util/ folder allows to extract certain live data from
!!   the output file. Run it with option -h for usage information.
!! - The script plot_live_data.sh allows to plot live data, e.g., the energy
!!   time traces. Run it with option -h for usage information.
!!
module live_data
  
#include "version.h"
  
  implicit none
  
  private
  public init_live_data, write_live_data, write_live_data_vacuum, finalize_live_data
  
  integer,           parameter :: LIVE_DATA_HANDLE = 43 !< File handle for live data file
  character(len=20), parameter :: LIVE_DATA_FILE   = 'macroscopic_vars.dat' !< Live data file
  
  
  
  contains
  
  
  
  !> Open file, write out headers and some parameters.
  subroutine init_live_data()
    
    use mod_parameters,    only: n_tor, n_plane, n_period, jorek_model, variable_names
    use phys_module,   only: produce_live_data, mode, mode_type, xpoint, xcase
    
    implicit none
    
    logical :: opened
    integer :: n, i
    
    if ( .not. produce_live_data ) return
    
    ! --- Check, that the file handle is not already in use.
    inquire(unit=LIVE_DATA_HANDLE, opened=opened)
    if ( opened ) then
      write(*,*) 'WARNING: LIVE DATA CANNOT BE PRODUCED AS FILE HANDLE IS ALREADY IN USE!'
      produce_live_data = .false.
      return
    end if
    
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='REPLACE', action='WRITE')
    
    ! --- Write some general information
    write(LIVE_DATA_HANDLE,*)  '@rcs_version: ', RCS_VERSION
    write(LIVE_DATA_HANDLE,'(A,I5)') '@jorek_model: ', jorek_model
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_tor: ', n_tor
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_plane: ', n_plane
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_period: ', n_period
    write(LIVE_DATA_HANDLE,'(A)') '@plottable: energies growth_rates times input_profiles axis current betas particlecontent thermalenergy heatingpower &
                                    particlesource diag_coil_curr integrated_energies bnd_fluxes dEdt helicity dissipative_terms work_terms             &
                                    area volume li3 energy_conservation'
    write(LIVE_DATA_HANDLE,'(A,15(A11,1X))') '@variable_names: ', variable_names
    
    ! --- Write file headers indicating what data is in the files.
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_times: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@times_xlabel: time step'
    write(LIVE_DATA_HANDLE,'(A)') '@times_ylabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@times_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@times: "step"     "time"'
    
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_energies: ', 2*(n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@energies_ylabel: normalized energy'
    write(LIVE_DATA_HANDLE,'(A)') '@energies_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@energies: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"E_{mag', mode(n), '}"'
    end do
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"E_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_growth_rates: ', 2*(n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_ylabel: normalized growth rate'
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@growth_rates: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{mag', mode(n), '}"'
    end do
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_axis: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@axis_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@axis_ylabel: Magnetic axis properties'
    write(LIVE_DATA_HANDLE,'(A)') '@axis_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@axis: %"time"           "R position"              "Z position"           "Psi on axis"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_integrated_energies: ', 5 
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_ylabel: Total integrated energies [J]'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies: %"time"           "Total energy"              "Magnetic"           "Kinetic parallel"    &
                                   "Kinetic perpendicular"                 "Thermal energy"     '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_bnd_fluxes: ', 4 
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_ylabel: Total boundary fluxes [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes: %"time"       "p vn"  "kinpar-flux"    "qn-par"    "qn-perp"   '
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_current: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@current_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@current_ylabel: plasma current [A]'
    write(LIVE_DATA_HANDLE,'(A)') '@current_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@current: %"time"       "Total"    "Inside LCFS"   "Outside LCFS" '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_helicity: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_ylabel: Total helicity [Wb^2]'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity: %"time"           "Helicity"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_area: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@area_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@area_ylabel: Total area [m^2]'
    write(LIVE_DATA_HANDLE,'(A)') '@area_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@area: %"time"           "Plasma cross sectional area"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_volume: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@volume_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@volume_ylabel: Total volume [m^3]'
    write(LIVE_DATA_HANDLE,'(A)') '@volume_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@volume: %"time"           "Plasma volume"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_li3: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@li3_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@li3_ylabel: li(3) '
    write(LIVE_DATA_HANDLE,'(A)') '@li3_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@li3: %"time"         "inside separatrix"  "All domain" '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_energy_conservation: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_ylabel: Total energy conservation'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation: %"time"       "-dEtotdt"     "Sum bnd fluxes + sources + dissipative terms" '
    write(LIVE_DATA_HANDLE,*)
 
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_dissipative_terms: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_ylabel: Dissipative powers [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms: %"time"         "Ohmic power"   "Parallel viscosity power" '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_work_terms: ', 2 
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_ylabel: Total work [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms: %"time"     "Magnetic = JxB~nabla p"   "Thermal = vpar*nabla p"  '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_dEdt: ', 5 
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_ylabel: -dEnergydt [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt: %"time"    "Etot"  "Wmagtot"  "thermaltot"   "kinperptot"  "kinpartot"      '
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_betas: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@betas_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@betas_ylabel: plasma beta'
    write(LIVE_DATA_HANDLE,'(A)') '@betas_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@betas: %"time"           "beta poloidal"       "beta toroidal"       "beta normalized"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_particlecontent: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_ylabel: particle content'
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent: %"time"           "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_thermalenergy: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_ylabel: thermal energy'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy: %"time"   "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_heatingpower: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_ylabel: heating power [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower: %"time"  "Total"     "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_particlesource: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_ylabel: particle source [10^20/m^3/s]'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource: %"time"  "Total"   "inside separatrix"   "outside separatrix"   '
    write(LIVE_DATA_HANDLE,*)
    
    ! --- Call the model-specific part of the init_live_data routine
    call init_live_data_model(LIVE_DATA_HANDLE) 
    
    close(LIVE_DATA_HANDLE)
    
  end subroutine init_live_data
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data(index)
    
    use mod_parameters,  only: n_tor
    use phys_module, only: xtime, energies, produce_live_data, R_axis_t, Z_axis_t, Psi_axis_t,     &
      current_t, beta_p_t, beta_t_t, beta_n_t, density_in_t, density_out_t, pressure_in_t,               &
      pressure_out_t, heat_src_in_t, heat_src_out_t, part_src_in_t, part_src_out_t, &
      E_tot_t, Helicity_tot_t, Kin_perp_tot_t, thermal_tot_t, kin_par_tot_t, ohmic_tot_t,      &
      Wmag_tot_t, Ip_tot_t, flux_pvn_t, flux_qpar_t, flux_qperp_t, flux_kinpar_t, dE_tot_dt, &
      dWmag_tot_dt, dthermal_tot_dt, dkinpar_tot_dt, dkinperp_tot_dt,  Magwork_tot_t,   &
      thmwork_tot_t, viscopar_dissip_tot_t, viscopar_flux_t, li3_t,      &
      li3_tot_t, part_src_tot_t, heat_src_tot_t, volume_t, area_t 


    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    integer :: i, j
    real*8  :: e1, e2, growth_rate, sum_fluxes_dissip
    
    if ( .not. produce_live_data ) return
    
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='OLD', position='APPEND', action='WRITE')
    
    ! --- Write data to the files.
    write(LIVE_DATA_HANDLE,'(A,I6,1X,ES17.9)') '@times:', index, xtime(index)
    write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@energies:', xtime(index)
    do j = 1, 2
      do i = 1, n_tor, 2
        write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') sum(energies(max(i-1,1):i,j,index))
      end do
    end do
    write(LIVE_DATA_HANDLE,*)
    if ( index > 1 ) then
      write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@growth_rates:', &
        (xtime(index)+xtime(index-1))/2.d0
      do j = 1, 2
        do i = 1, n_tor, 2
          e1 = sum(energies(max(i-1,1):i,j,index))
          e2 = sum(energies(max(i-1,1):i,j,index-1))
          if ( (e1 .NE. 0.) .and. (e2 .NE. 0.) ) then
             growth_rate = 0.5d0 * ( log(e1) - log(e2) ) / (xtime(index)-xtime(index-1))
          else
             growth_rate = 0.d0
          endif
          write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') growth_rate
        end do
      end do
    end if
    write(LIVE_DATA_HANDLE,*)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@axis: ', xtime(index), R_axis_t(index), Z_axis_t(index), Psi_axis_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@current: ', xtime(index), Ip_tot_t(index), current_t(index), Ip_tot_t(index)-current_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@betas: ', xtime(index), beta_p_t(index), beta_t_t(index), beta_n_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlecontent: ', xtime(index), density_in_t(index), density_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@thermalenergy: ', xtime(index), pressure_in_t(index), pressure_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@heatingpower: ', xtime(index), heat_src_tot_t(index), heat_src_in_t(index), heat_src_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlesource: ', xtime(index), part_src_tot_t(index), part_src_in_t(index), part_src_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@integrated_energies: ', xtime(index), E_tot_t(index), Wmag_tot_t(index), &
                                                     kin_par_tot_t(index),  kin_perp_tot_t(index),  thermal_tot_t(index)  
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@helicity: ', xtime(index), helicity_tot_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@area: ', xtime(index), area_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@volume: ', xtime(index), volume_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@li3: ', xtime(index), li3_t(index), li3_tot_t(index)

    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dissipative_terms: ', xtime(index), ohmic_tot_t(index), viscopar_dissip_tot_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@work_terms: ', xtime(index), Magwork_tot_t(index), thmwork_tot_t(index)
    write(LIVE_DATA_HANDLE,'(A,7ES17.9)') '@bnd_fluxes: ', xtime(index), flux_Pvn_t(index), flux_kinpar_t(index), &
                                           flux_qpar_t(index), flux_qperp_t(index)

   if(index>1) then
       write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@dEdt: ', xtime(index-1), -dE_tot_dt(index-1), -dWmag_tot_dt(index-1), &
                                            -dthermal_tot_dt(index-1),-dkinperp_tot_dt(index-1),-dkinpar_tot_dt(index-1)
   
       sum_fluxes_dissip = flux_Pvn_t(index-1)  + flux_kinpar_t(index-1) + flux_qpar_t(index-1) + flux_qperp_t(index-1) &
                         + ohmic_tot_t(index-1) + viscopar_dissip_tot_t(index-1) - heat_src_tot_t(index-1)
 
      write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@energy_conservation: ', xtime(index-1), -dE_tot_dt(index-1), sum_fluxes_dissip 

    else
      write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@dEdt: ', xtime(index), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@energy_conservation: ', xtime(index), 0.d0, 0.d0 
    endif
 
    close(LIVE_DATA_HANDLE)
    
  end subroutine write_live_data
  
  
  
  !> Close file.
  subroutine finalize_live_data()
    
    use phys_module, only: produce_live_data
    
    implicit none
    
    if ( .not. produce_live_data ) return
    
    ! -nothing to be done currently-
    
  end subroutine finalize_live_data
  
  
  
  subroutine write_live_data_vacuum(index, diag_coil_curr)
    
    use phys_module, only: xtime
    
    integer,             intent(in) :: index
    real*8, allocatable, intent(in) :: diag_coil_curr(:,:)
    
    logical, save :: header_written = .false.
    
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='OLD', position='APPEND', action='WRITE')
    
    if ( allocated(diag_coil_curr) ) then
      if ( .not. header_written ) then
        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_diag_coil_curr: ', size(diag_coil_curr,2)
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_ylabel: Diagnostic coil current'
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_logy: 0'
        write(LIVE_DATA_HANDLE,*)
        header_written = .true.
      end if
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@diag_coil_curr: ', xtime(index), diag_coil_curr(index,:)
    end if
    
    close(LIVE_DATA_HANDLE)
    
  end subroutine write_live_data_vacuum
  
end module live_data

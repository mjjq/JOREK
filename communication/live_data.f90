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
 
    write(*,*) "CCC111", produce_live_data   
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
    write(LIVE_DATA_HANDLE,'(A)') '@plottable: energies growth_rates times input_profiles axis current betas particlecontent thermalenergy heatingpower particlesource diag_coil_curr'
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
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_current: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@current_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@current_ylabel: plasma current'
    write(LIVE_DATA_HANDLE,'(A)') '@current_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@current: %"time"           "Current"'
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
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_heatingpower: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_ylabel: heating power'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower: %"time"   "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_particlesource: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_ylabel: particle source'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource: %"time"   "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(*,*) "CCC222"

    ! --- Call the model-specific part of the init_live_data routine
    call init_live_data_model(LIVE_DATA_HANDLE) 
    
      write(*,*) "CCC333"

    close(LIVE_DATA_HANDLE)
    
  end subroutine init_live_data
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data(index)
    
    use mod_parameters,  only: n_tor
    use phys_module, only: xtime, energies, produce_live_data, R_axis_t, Z_axis_t, Psi_axis_t,     &
      current_t, beta_p_t, beta_t_t, beta_n_t, density_in_t, density_out_t, pressure_in_t,               &
      pressure_out_t, heat_src_in_t, heat_src_out_t, part_src_in_t, part_src_out_t
    
    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    integer :: i, j
    real*8  :: e1, e2, growth_rate
    
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
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@current: ', xtime(index), current_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@betas: ', xtime(index), beta_p_t(index), beta_t_t(index), beta_n_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlecontent: ', xtime(index), density_in_t(index), density_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@thermalenergy: ', xtime(index), pressure_in_t(index), pressure_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@heatingpower: ', xtime(index), heat_src_in_t(index), heat_src_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlesource: ', xtime(index), part_src_in_t(index), part_src_out_t(index)
    
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

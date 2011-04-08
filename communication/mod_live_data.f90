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
  
  implicit none
  
  private
  public init_live_data, write_live_data, finalize_live_data
  
  integer,           parameter :: LIVE_DATA_HANDLE = 43 !< File handle for live data file
  character(len=20), parameter :: LIVE_DATA_FILE   = 'macroscopic_vars.dat' !< Live data file
  
  
  
  contains
  
  
  
  !> Open file, write out headers and some parameters.
  subroutine init_live_data()
    
    use parameters,  only: n_tor, n_plane, n_period, jorek_model, variable_names
    use phys_module, only: produce_live_data, mode, mode_type
    
    implicit none
    
    logical :: opened
    integer :: n
    
    if ( .not. produce_live_data ) return
    
    ! --- Check, that the file handle is not already in use.
    inquire(unit=LIVE_DATA_HANDLE, opened=opened)
    if ( opened ) then
      write(*,*) 'WARNING: LIVE DATA CANNOT BE PRODUCED AS FILE HANDLE IS ALREADY IN USE!'
      produce_live_data = .false.
      return
    end if
    
    ! --- Open the data files.
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='REPLACE', action='WRITE')
    
    ! --- Write some general information
    write(LIVE_DATA_HANDLE,'(A,I5)') '@jorek_model: ', jorek_model
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_tor: ', n_tor
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_plane: ', n_plane
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_period: ', n_period
    write(LIVE_DATA_HANDLE,'(A)') '@plottable: energies growth_rates times'
    write(LIVE_DATA_HANDLE,'(A,15(A11,1X))') '@variable_names: ', variable_names
    
    ! --- Write file headers indicating what data is in the files.
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_times: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@times_xlabel: time step'
    write(LIVE_DATA_HANDLE,'(A)') '@times_ylabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@times: "step"     "time"'
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_energies: ', 2*n_tor
    write(LIVE_DATA_HANDLE,'(A)') '@energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@energies_ylabel: normalized energy'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@energies: "time"           '
    do n = 1, n_tor
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,",",A3,A2,1x)',advance='no') '"E_{mag', mode(n), &
        mode_type(n), '}"'
    end do
    do n = 1, n_tor
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,",",A3,A2,1x)',advance='no') '"E_{kin', mode(n), &
        mode_type(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_growth_rates: ', 2*n_tor
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_ylabel: normalized growth rate'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@growth_rates: "time"           '
    do n = 1, n_tor
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,",",A3,A2,1x)',advance='no') '"G_{mag', mode(n), &
        mode_type(n), '}"'
    end do
    do n = 1, n_tor
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,",",A3,A2,1x)',advance='no') '"G_{kin', mode(n), &
        mode_type(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)
    
  end subroutine init_live_data
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data(index)
    
    use parameters,  only: n_tor
    use phys_module, only: xtime, energies, produce_live_data
    
    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    real*8 :: growth_rates(1:n_tor,1:2)
    
    if ( .not. produce_live_data ) return
    
    ! --- Write data to the files.
    write(LIVE_DATA_HANDLE,'(A,I6,1X,ES17.9)') '@times:', index, xtime(index)
    write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@energies:', xtime(index), energies(1:n_tor,1:2,index)
    if ( index > 1 ) then
      where ( (energies(:,:,index)>0.d0) .and. (energies(:,:,index-1)>0.d0) )
        growth_rates =                                                                         &
          0.5d0 * ( LOG(energies(1:n_tor,1:2,index)) - LOG(energies(1:n_tor,1:2,index-1)) )    &
          / (xtime(index)-xtime(index-1))
      elsewhere
        growth_rates = 0.d0
      end where
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@growth_rates:', &
        (xtime(index)+xtime(index-1))/2.d0, growth_rates
    end if
    
  end subroutine write_live_data
  
  
  
  !> Close file.
  subroutine finalize_live_data()
    
    use phys_module, only: produce_live_data
    
    implicit none
    
    if ( .not. produce_live_data ) return
    
    close(LIVE_DATA_HANDLE)
    
  end subroutine finalize_live_data

end module live_data

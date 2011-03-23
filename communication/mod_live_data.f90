!> The module contains routines to write certain data to text files while the code is running.
!!
!! The input parameter phys_module::produce_live_data allows to switch the functionality of
!! this module on or off. If the parameter is true, the following files are created:
!! <ul><li>
!!   plot_energies.gnuplot: Allows to plot the energy time evolution of the n=0,
!!   and the highest mode number using \n <tt>gnuplot plot_energies.gnuplot</tt>
!! </li><li>
!!   energies.dat: Magnetic and kinetic energies versus time
!! </li><li>
!!   growth_rates.dat: Growth rates versus time
!! </li><li>
!!   times.dat: JOREK time versus time step index
!! </li></ul>
!!
!! The data can be plotted while or after JOREK is running to monitor a simulation.
!! E.g., start gnuplot and type: \n
!!     <tt> set log y; plot 'energies.dat' using 1:3 with linespoints </tt>
!!
module live_data
  
  implicit none
  
  private
  public init_live_data, write_live_data, finalize_live_data
  
  integer, parameter :: TIMES_FILE    = 43 !< File handle for 'times.dat'
  integer, parameter :: ENERGIES_FILE = 44 !< File handle for 'energies.dat'
  integer, parameter :: GROWTH_FILE   = 45 !< File handle for 'growth_rates.dat'
  
  
  
  contains
  
  
  
  !> Open files, the data is written to.
  subroutine init_live_data()
    
    use parameters,  only: n_tor
    use phys_module, only: produce_live_data
    
    implicit none
    
    logical :: opened1, opened2, opened3
    
    if ( .not. produce_live_data ) return
    
    ! --- Check, that file handles are not already in use.
    inquire(unit=TIMES_FILE, opened=opened1)
    inquire(unit=ENERGIES_FILE, opened=opened2)
    inquire(unit=GROWTH_FILE, opened=opened3)
    if ( opened1 .or. opened2 .or. opened3 ) then
      write(*,*) 'WARNING: LIVE DATA CANNOT BE PRODUCED AS A FILE HANDLE IS ALREADY IN USE!'
      write(*,*) opened1, opened2, opened3
      produce_live_data = .false.
      return
    end if
    
    ! --- Open the data files.
    open(TIMES_FILE,    file='times.dat',        status='REPLACE', action='WRITE')
    open(ENERGIES_FILE, file='energies.dat',     status='REPLACE', action='WRITE')
    open(GROWTH_FILE,   file='growth_rates.dat', status='REPLACE', action='WRITE')
    
  end subroutine init_live_data
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data(index)
    
    use phys_module, only: xtime, energies, produce_live_data
    
    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    if ( .not. produce_live_data ) return
    
    ! --- Write data to the files.
    write(TIMES_FILE,'(I6,1X,ES15.7)') index, xtime(index)
    write(ENERGIES_FILE,'(999ES15.7)') xtime(index), energies(:,1:2,index)
    if ( index > 1 ) then 
      write(GROWTH_FILE,  '(999ES15.7)') (xtime(index)+xtime(index-1))/2.d0,             &
        0.5d0 * ( LOG(energies(:,1:2,index)) - LOG(energies(:,1:2,index-1)) )            &
        / (xtime(index)-xtime(index-1))
    end if
    
  end subroutine write_live_data
  
  
  
  !> Close files, the data is written to.
  subroutine finalize_live_data()
    
    use phys_module, only: produce_live_data
    
    implicit none
    
    if ( .not. produce_live_data ) return
    
    ! --- Close the data files.
    close(TIMES_FILE)
    close(ENERGIES_FILE)
    close(GROWTH_FILE)
    
  end subroutine finalize_live_data

end module live_data

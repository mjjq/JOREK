!> The module contains routines to write certain data to text files while the code is running.
!!
!! The input parameter phys_module::produce_live_data allows to switch the functionality of
!! this module on or off. If the parameter is true, the following files are created:
!! - energies.dat: Magnetic and kinetic energies versus time
!! - growth_rates.dat: Growth rates versus time
!! - times.dat: JOREK time versus time step index
!!
!! The data can be plotted while or after JOREK is running to monitor a simulation.
!! E.g., start gnuplot and type: \n
!!     <tt> set log y; set nokey; plot for [i=2:7] 'energies.dat' using 1:i with linespoints </tt>
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
    use phys_module, only: produce_live_data, mode, mode_type
    
    implicit none
    
    logical :: opened1, opened2, opened3
    integer :: n
    
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
    
    ! --- Write file headers indicating what data is in the files.
    write(TIMES_FILE,'(A)') '# step       time'

    write(ENERGIES_FILE,'(A)',advance='no') '#      time      '
    do n = 1, n_tor
      write(ENERGIES_FILE,'(A5,",",I2,",",A3,3x)',advance='no') 'E_mag', mode(n), mode_type(n)
    end do
    do n = 1, n_tor
      write(ENERGIES_FILE,'(A5,",",I2,",",A3,3x)',advance='no') 'E_kin', mode(n), mode_type(n)
    end do
    write(ENERGIES_FILE,*)

    write(GROWTH_FILE,'(A)',advance='no') '#      time      '
    do n = 1, n_tor
      write(GROWTH_FILE,'(A5,",",I2,",",A3,3x)',advance='no') 'G_mag', mode(n), mode_type(n)
    end do
    do n = 1, n_tor
      write(GROWTH_FILE,'(A5,",",I2,",",A3,3x)',advance='no') 'G_kin', mode(n), mode_type(n)
    end do
    write(GROWTH_FILE,*)
    
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
    write(TIMES_FILE,'(I6,1X,ES15.7)') index, xtime(index)
    write(ENERGIES_FILE,'(999ES15.7)') xtime(index), energies(1:n_tor,1:2,index)
    if ( index > 1 ) then
      where ( (energies(:,:,index)>0.d0) .and. (energies(:,:,index-1)>0.d0) )
        growth_rates =                                                                         &
          0.5d0 * ( LOG(energies(1:n_tor,1:2,index)) - LOG(energies(1:n_tor,1:2,index-1)) )    &
          / (xtime(index)-xtime(index-1))
      elsewhere
        growth_rates = 0.d0
      end where
      write(GROWTH_FILE,  '(999ES15.7)') (xtime(index)+xtime(index-1))/2.d0, growth_rates
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

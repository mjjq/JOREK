!> The module contains routines to write certain data to text files during the code run
!! (e.g., energies and growth rates).
module live_data
  
  implicit none
  
  private
  public init_live_data, write_live_data, finalize_live_data
  
  integer, parameter :: TIMES_FILE    = 43 !< File handler for 'times.dat'
  integer, parameter :: ENERGIES_FILE = 44 !< File handler for 'energies.dat'
  integer, parameter :: GROWTH_FILE   = 45 !< File handler for 'growth_rates.dat'
  
  
  
  contains
  
  
  
  !> Open files, the data is written to.
  subroutine init_live_data()
    
    implicit none
    
    open(TIMES_FILE,    file='times.dat',        status='REPLACE', action='WRITE')
    open(ENERGIES_FILE, file='energies.dat',     status='REPLACE', action='WRITE')
    open(GROWTH_FILE,   file='growth_rates.dat', status='REPLACE', action='WRITE')
    
  end subroutine init_live_data
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data(index)
    
    use phys_module, only: xtime, energies
    
    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    write(TIMES_FILE,'(I6,1X,ES13.5)') index, xtime(index)
    write(ENERGIES_FILE,'(999ES13.5)') xtime(index), energies(:,1:2,index)
    if ( index > 1 ) then 
      write(GROWTH_FILE,  '(999ES13.5)') (xtime(index)+xtime(index-1))/2.d0,             &
        0.5d0 * ( LOG(energies(:,1:2,index)) - LOG(energies(:,1:2,index-1)) )            &
        / (xtime(index)-xtime(index-1))
    end if
    
  end subroutine write_live_data
  
  
  
  !> Close files, the data is written to.
  subroutine finalize_live_data()
    
    implicit none
    
    close(TIMES_FILE)
    close(ENERGIES_FILE)
    close(GROWTH_FILE)
    
  end subroutine finalize_live_data

end module live_data

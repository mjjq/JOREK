module matio_module
  use HDF5          
  Implicit None

  logical :: fileexist=.false.

  private
  public :: read_matrix_h5, save_mat_h5, timestamp, slurmid

!  interface timestamp
!    module procedure timestamp
!  end interface timestamp

  contains 


   subroutine read_matrix_h5(fname,n,nnz,rowptr,colptr,val,rhs)
      integer :: i,stat, ierr
      integer(HID_T) fid
!      real(kind=c_double),  dimension(:), allocatable, target, intent(OUT)    :: rhs        

      character(len=255),intent(in) :: fname
      integer,intent(out) :: nnz,n
      integer :: indexing

      integer, dimension(:), pointer :: rowptr
      integer, dimension(:), pointer :: colptr
      real(kind=c_double), dimension(:), pointer    :: val
      real(kind=c_double), dimension(:), pointer    :: rhs   

      write(*,*) "Reading ",trim(fname)

      call HDF5_open(trim(fname),fid,ierr)
      if (ierr/=0) return
      call HDF5_integer_reading(fid,n,"n")  
      call HDF5_integer_reading(fid,nnz,"nnz")        
      call HDF5_integer_reading(fid,indexing,"indexing")              
     
      allocate(colptr(nnz),rowptr(nnz))      
      allocate(val(nnz))
      allocate(rhs(n))

      call HDF5_array1D_reading_int(fid,rowptr,"irn")
      call HDF5_array1D_reading_int(fid,colptr,"jcn")      
      call HDF5_array1D_reading(fid,val,"val")            
      call HDF5_array1D_reading(fid,rhs,"rhs")
      call HDF5_close(fid)

! change to C indexing      
      if (indexing==1) then
              do i=1,nnz
                rowptr(i)=rowptr(i)-1
                colptr(i)=colptr(i)-1                
              enddo
      endif


      return
  end subroutine read_matrix_h5  
! ########################### save global matrix to HD5 file ############################
  subroutine save_mat_h5(rank,n,nnz,irn,jcn,val,rhs)
    use hdf5_io_module
    
    integer :: rank,n,nnz,ierr
    integer(HID_T) fid
    CHARACTER(LEN=10)               :: fname
    integer, dimension(:), pointer       :: irn,jcn
    real(kind=c_double)    , dimension(:), pointer       :: val
    real(kind=c_double), dimension(:), pointer, optional       :: rhs

    write(fname,'(A5,I2.2,A3)') "matA_",rank,".h5"
    

    call HDF5_create(filename=fname,file_id=fid,ierr=ierr)
    call HDF5_integer_saving(fid,n,'n')
    call HDF5_integer_saving(fid,nnz,'nnz')
    call HDF5_integer_saving(fid,1,'indexing')  
    call HDF5_array1D_saving_int(fid,irn,nnz,'irn')
    call HDF5_array1D_saving_int(fid,jcn,nnz,'jcn')   
    call HDF5_array1D_saving(fid,val,nnz,'val')
    if (present(rhs)) call HDF5_array1D_saving(fid,rhs,n,'rhs')
    call HDF5_close(fid)
    
  end subroutine save_mat_h5  
!#############################################################################################
   subroutine save_solution_h5(fname,n,x)
      integer(HID_T) fid
      integer :: i,stat, ierr
      integer, intent(in) :: n
      character(len=255),intent(in) :: fname
      real*8, dimension(:), intent(in) :: x

      write(*,*) "Saving ",trim(fname)

      call HDF5_create(trim(fname),fid,ierr)
      call HDF5_integer_saving(fid,n,"n")  
      call HDF5_array1D_saving(fid,x,n,"x")            
      call HDF5_close(fid)        

      return
  end subroutine save_solution_h5  
!#############################################################################################
  subroutine HDF5_open(filename,file_id,ierr)
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr

    integer :: ierr_HDF5

    call H5open_f(ierr_HDF5)
    call H5Fopen_f(trim(filename)//char(0),H5F_ACC_RDONLY_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_open

  !----------------------------------------
  ! create HDF5 file 
  !----------------------------------------
  subroutine HDF5_create(filename,file_id,ierr)
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr

    integer :: ierr_HDF5

    call H5open_f(ierr_HDF5)
    call H5Fcreate_f(trim(filename)//char(0),H5F_ACC_TRUNC_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_create
  !----------------------------------------
  ! close HDF5 file 
  !----------------------------------------
  subroutine HDF5_close(file_id)
    integer(HID_T), intent(in) :: file_id   ! file identifier
    integer :: error   ! error flag

    call H5Fclose_f(file_id,error)
  end subroutine HDF5_close  
  !----------------------------------------
  ! HDF5 file: read integer
  !----------------------------------------
  subroutine HDF5_integer_reading(file_id,int,dsetname)
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    integer         , intent(out) :: int
    character(LEN=*), intent(in)  :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type

    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_reading
  !----------------------------------------
  ! HDF5 file: read 1D integer array
  !----------------------------------------
  subroutine HDF5_array1D_reading_int(file_id,array1D,dsetname)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(out), dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_int  

  !---------------------------------------- 
  ! HDF5 reading of 1D array of real*8
  !----------------------------------------
  subroutine HDF5_array1D_reading(file_id,array1D,dsetname)
    integer(HID_T), intent(in)   :: file_id   ! file identifier
    real*8        , dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading 

  !---------------------------------------- 
  ! HDF5 saving for a integer
  !----------------------------------------
  subroutine HDF5_integer_saving(file_id,int,dsetname)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(in) :: int
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier

    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,dataspace,dataset,error)
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_saving


  !---------------------------------------- 
  ! HDF5 saving for a real double
  !----------------------------------------
  subroutine HDF5_real_saving(file_id,rd,dsetname)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    real*8          , intent(in) :: rd
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier

    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)
    call H5Dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dataspace,dataset,error)
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_saving


  !---------------------------------------- 
  ! HDF5 saving for a 1D array of integer
  !----------------------------------------
  subroutine HDF5_array1D_saving_int(file_id,array1D,dim1,dsetname)
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    integer, dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,filespace,dataset,error)
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_int  

!-------------------------------------------- 
  ! HDF5 saving of 1D array of real*8
  !--------------------------------------------
  subroutine HDF5_array1D_saving(file_id,array1D,dim1,dsetname)
    integer(HID_T)      , intent(in) :: file_id   ! file identifier
    real*8, dimension(:), intent(in) :: array1D
    integer             , intent(in) :: dim1
    character(LEN=*)    , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE,filespace,dataset,error)
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving  
!###########################################################################################      
  subroutine slurmid(rank)
    character(len=12) :: envname="SLURM_PROCID"
    character(len=4) :: val
    integer :: rank

    call get_environment_variable (envname, val)
    open(unit = 101, file = 'procid.out', status='REPLACE', action='WRITE')
    write(101,*) val, rank
    close(101)

  end subroutine slurmid

  subroutine timestamp(msg,id)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: id

    integer,dimension(8) :: values
    character(len=14) :: fname
    real :: t 
    
    if (present(id)) then
      write (fname, "(A8,(I0.2),A4)") 'timeline', id, '.out'
    else
      write (fname, "(A12)") 'timeline.out'
    endif

    if (.not.fileexist) open(unit = 100, file = trim(fname), status='REPLACE', action='WRITE')
    fileexist = .true.

    call date_and_time(VALUES=values)
    t = values(5)*3600000 + values(6)*60000 + values(7)*1000 + values(8)
    open(unit = 100, file = trim(fname), status='OLD', position="append", action='WRITE')
    write(100,*) t, trim(msg) 
    close(100)

  end subroutine timestamp


end module matio_module

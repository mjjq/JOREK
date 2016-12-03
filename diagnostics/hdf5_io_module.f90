!---------------------------------------------
! file : HDF5_io.f90
! date : 25/01/2006
!  array saving and reading in HDF5 format
!
! Adapted from GYSELA routines
!---------------------------------------------
module hdf5_io_module
  implicit none

#ifdef USE_HDF5  

  !******************************
  contains
  !******************************

  !---------------------------------------- 
  ! create HDF5 file 
  !----------------------------------------
  subroutine HDF5_create(filename,file_id,ierr)
    use HDF5
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr

    integer :: ierr_HDF5

    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)

    !*** Create a new file using default properties ***
    call H5Fcreate_f(trim(filename)//char(0), &
      H5F_ACC_TRUNC_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_create


  !---------------------------------------- 
  ! open HDF5 file 
  !----------------------------------------
  subroutine HDF5_open(filename,file_id,ierr)
    use HDF5
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr

    integer :: ierr_HDF5

    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)

    !*** open the HDF5 file ***
    call H5Fopen_f(trim(filename)//char(0), &
      H5F_ACC_RDONLY_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_open


  !---------------------------------------- 
  ! close HDF5 file 
  !----------------------------------------
  subroutine HDF5_close(file_id)
    use HDF5
    integer(HID_T), intent(in) :: file_id   ! file identifier

    integer :: error   ! error flag

    call H5Fclose_f(file_id,error)
  end subroutine HDF5_close


  !*************************************************
  !  HDF5 WRITING
  !*************************************************
  !---------------------------------------- 
  ! HDF5 saving for a character string
  !----------------------------------------
  subroutine HDF5_char_saving(file_id,charvar,dsetname)
    use HDF5
    use H5LT
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    character(LEN=*), intent(in) :: charvar
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: ierr_HDF5  ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: type_id  ! dataspace identifier

   !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    !*** Create character dataset ***
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,ierr_HDF5)
    CALL h5tset_size_f(type_id,int(len(charvar),SIZE_T),ierr_HDF5)
    CALL h5screate_f(H5S_SCALAR_F,dataspace,ierr_HDF5)
    CALL h5dcreate_f(file_id,dsetname,type_id,dataspace,dataset,ierr_HDF5)
    CALL h5dwrite_f(dataset,type_id,charvar,dim,ierr_HDF5)
    !*** Closing ***
    call H5Dclose_f(dataset,ierr_HDF5)
    call H5Sclose_f(dataspace,ierr_HDF5)
  end subroutine HDF5_char_saving


  !---------------------------------------- 
  ! HDF5 saving for a character 1D array
  !----------------------------------------
  subroutine HDF5_array1D_saving_char(file_id,array1D,dim1,dsetname)
    use HDF5
    use H5LT
    integer(HID_T)                , intent(in) :: file_id   ! file identifier
    character(LEN=*), dimension(:), intent(in) :: array1D
    integer                       , intent(in) :: dim1
    character(LEN=*)              , intent(in) :: dsetname  ! dataset name

    integer             :: ierr_HDF5  ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,ierr_HDF5)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_CHARACTER,&
      dataspace,dataset,ierr_HDF5)

    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    call H5Dwrite_f(dataset,H5T_NATIVE_CHARACTER,array1D,dim,ierr_HDF5)

    !*** Closing ***
    call H5Sclose_f(dataspace,ierr_HDF5)
    call H5Dclose_f(dataset,ierr_HDF5)
  end subroutine HDF5_array1D_saving_char


  !---------------------------------------- 
  ! HDF5 saving for an integer 
  !----------------------------------------
  subroutine HDF5_integer_saving(file_id,int,dsetname)
    use HDF5
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(in) :: int
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,trim(dsetname), &
      H5T_NATIVE_INTEGER,dataspace,dataset,error)

    !*** Write the integer data to the dataset ***
    !***  using default transfer properties    ***
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)

    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_saving


  !---------------------------------------- 
  ! HDF5 saving for a real double
  !----------------------------------------
  subroutine HDF5_real_saving(file_id,rd,dsetname)
    use HDF5
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    real*8          , intent(in) :: rd
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,dsetname, &
      H5T_NATIVE_DOUBLE,dataspace,dataset,error)

    !*** Write the integer data to the dataset ***
    !***  using default transfer properties    ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)

    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_saving


  !---------------------------------------- 
  ! HDF5 saving for a 1D array of integer
  !----------------------------------------
  subroutine HDF5_array1D_saving_int(file_id,array1D,dim1,dsetname)
    use HDF5
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
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      dataspace,dataset,error)

    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)

    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_int

  !---------------------------------------- 
  ! HDF5 saving for a 2D array of integer
  !----------------------------------------
  subroutine HDF5_array2D_saving_int(file_id,array2D,dim1,dim2,dsetname)
    use HDF5
    integer(HID_T)           , intent(in) :: file_id   ! file identifier
    integer, dimension(:,:)  , intent(in) :: array2D
    integer                  , intent(in) :: dim1, dim2
    character(LEN=*)         , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      dataspace,dataset,error)

    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error)

    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_saving_int


  !---------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*4
  !----------------------------------------
  subroutine HDF5_array1D_saving_r4(file_id,array1D,dim1,dsetname)
    use HDF5
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    real(4), dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_REAL,&
      dataspace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_REAL,array1D,dim,error)

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_r4


  !-------------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*8
  !--------------------------------------------
  subroutine HDF5_array1D_saving(file_id,array1D,dim1,dsetname)
    use HDF5
    integer(HID_T)      , intent(in) :: file_id   ! file identifier
    real*8, dimension(:), intent(in) :: array1D
    integer             , intent(in) :: dim1
    character(LEN=*)    , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Creates a new property for gzip dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE,&
      dataspace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 2D array
  !----------------------------------------
  subroutine HDF5_array2D_saving(file_id,array2D,dim1,dim2,dsetname)
    use HDF5
    integer(HID_T)        , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:), intent(in) :: array2D
    integer               , intent(in) :: dim1, dim2
    character(LEN=*)      , intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer              :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(2)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 3D array
  !----------------------------------------
  subroutine HDF5_array3D_saving(file_id,array3D, &
    dim1,dim2,dim3,dsetname)
    use HDF5
    integer(HID_T)          , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:,:), intent(in) :: array3D
    integer                 , intent(in) :: dim1, dim2, dim3
    character(LEN=*)        , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    rank   = 3
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 4D array
  !----------------------------------------
  subroutine HDF5_array4D_saving(file_id,array4d, &
    dim1,dim2,dim3,dim4,dsetname)
    use HDF5
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:,:,:), intent(in) :: array4d
    integer                   , intent(in) :: dim1, dim2, dim3, dim4
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier
 
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    rank   = 4
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error)

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array4D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 5D array
  !----------------------------------------
  subroutine HDF5_array5D_saving(file_id,array5d, &
    dim1,dim2,dim3,dim4,dim5,dsetname)
    use HDF5
    integer(HID_T)              , intent(in) :: file_id  ! file identifier
    real*8, dimension(:,:,:,:,:), intent(in) :: array5d
    integer                     , intent(in) :: dim1, dim2
    integer                     , intent(in) :: dim3, dim4, dim5
    character(LEN=*)            , intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer             :: cmpr       ! compression level
    integer(HSIZE_T), &
      dimension(5)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    dim(5) = dim5
    rank   = 5
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)
    call H5Pset_chunk_f(property,rank,dim,error)
    cmpr = 6
    call H5Pset_deflate_f(property,cmpr,error)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      dataspace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties         ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_saving


  !************************************************
  !  HDF5 READING
  !************************************************
  !---------------------------------------- 
  ! HDF5 reading for a character
  !----------------------------------------
  subroutine HDF5_char_reading(file_id,charvar,dsetname)
    use HDF5
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    character       , intent(out) :: charvar
    character(LEN=*), intent(in)  :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_CHARACTER,charvar,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_char_reading


  !---------------------------------------- 
  ! HDF5 reading for a character array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading_char(file_id,array1D,dsetname)
    use HDF5
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    character       , intent(out), dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_CHARACTER,array1D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_char


  !---------------------------------------- 
  ! HDF5 reading for an integer 
  !----------------------------------------
  subroutine HDF5_integer_reading(file_id,int,dsetname)
    use HDF5
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

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_reading


  !---------------------------------------- 
  ! HDF5 reading for an integer array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading_int(file_id,array1D,dsetname)
    use HDF5
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(out), dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_int


  !---------------------------------------- 
  ! HDF5 reading for a real double
  !----------------------------------------
  subroutine HDF5_real_reading(file_id,rd,dsetname)
    use HDF5
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    real*8          , intent(out) :: rd
    character(LEN=*), intent(in)  :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_reading


  !---------------------------------------- 
  ! HDF5 reading for an array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading(file_id,array1D,dsetname)
    use HDF5
    integer(HID_T), intent(in)   :: file_id   ! file identifier
    real*8        , dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading


  !---------------------------------------- 
  ! HDF5 reading for an integer array 2D 
  !----------------------------------------
  subroutine HDF5_array2D_reading_int(file_id,array2D,dsetname)
    use HDF5
    integer(HID_T)  , intent(in)     :: file_id   ! file identifier
    integer         , dimension(:,:) :: array2D
    character(LEN=*), intent(in)     :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file opennning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array2D,1)
    dim(2) = size(array2D,2)
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading_int


  !---------------------------------------- 
  ! HDF5 reading for an array 2D
  !----------------------------------------
  subroutine HDF5_array2D_reading(file_id,array2D,dsetname)
    use HDF5
    integer(HID_T)  , intent(in)     :: file_id   ! file identifier
    real*8          , dimension(:,:) :: array2D
    character(LEN=*), intent(in)     :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file opennning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array2D,1)
    dim(2) = size(array2D,2)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading


  !---------------------------------------- 
  ! HDF5 reading for an array 3D
  !----------------------------------------
  subroutine HDF5_array3D_reading(file_id,array3D,dsetname,&
       in1, in2, in3)
    use HDF5
    integer(HID_T)     , intent(in)       :: file_id   ! file identifier
    real*8             , dimension(:,:,:) :: array3D
    character(LEN=*)   , intent(in)       :: dsetname  ! dataset name
    integer, intent(in), optional  :: in1, in2, in3
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    if (present(in1) .and. present(in2) .and. present(in3)) then
       dim(1) = in1
       dim(2) = in2
       dim(3) = in3
    else
       dim(1) = size(array3D,1)
       dim(2) = size(array3D,2)
       dim(3) = size(array3D,3)
    end if
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_reading


  !---------------------------------------- 
  ! HDF5 reading for an array 4D
  !----------------------------------------
  subroutine HDF5_array4D_reading(file_id,array4D,dsetname,ierr)
    use HDF5
    integer(HID_T), intent(in)     :: file_id   ! file identifier
    real*8, dimension(:,:,:,:)     :: array4D
    character(LEN=*), intent(in)   :: dsetname  ! dataset name
    integer, optional, intent(out) :: ierr

    integer             :: ierr_HDF5      ! error flag
    integer             :: rank           ! dataset rank
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,ierr_HDF5)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array4D,1)
    dim(2) = size(array4D,2)
    dim(3) = size(array4D,3)
    dim(4) = size(array4D,4)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,ierr_HDF5)

    !*** Closing ***
    call H5Dclose_f(dataset,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_array4D_reading


  !---------------------------------------- 
  ! HDF5 reading for an array 5D
  !----------------------------------------
  subroutine HDF5_array5D_reading(file_id,array5D,dsetname)
    use HDF5
    integer(HID_T), intent(in)   :: file_id   
    real*8, dimension(:,:,:,:,:) :: array5D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
           dimension(5) :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file openning ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim(1) = size(array5D,1)
    dim(2) = size(array5D,2)
    dim(3) = size(array5D,3)
    dim(4) = size(array5D,4)
    dim(5) = size(array5D,5)
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_reading

#endif
end module hdf5_io_module

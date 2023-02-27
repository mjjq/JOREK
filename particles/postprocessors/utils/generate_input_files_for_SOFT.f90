!> Generate magnetic equilibrium and distribution function
!> files compatible with the SOFT code inputs from the 
!> JOREK MHD solutions and distribution functions.
program generate_input_files_for_SOFT
use mod_mpi_tools,   only: init_mpi_threads,finalize_mpi_threads
use particle_tracer

implicit none

!> Variables --------------------------------------------------------------------------------------------
type(event)       :: field_reader
integer           :: my_id,n_cpus,ierr
integer           :: n_vec,n_R,n_Z,n_R_loc,n_Z_loc
real*8,dimension(:),allocatable     :: R_mesh,Z_mesh
real*8,dimension(:,:,:),allocatable :: magnetic_field
character(len=17) :: fields_filename
character(len=20) :: magnetic_field_filename
character(len=20) :: pdf_filename
!> Variables definitions --------------------------------------------------------------------------------
n_vec = 3 !< number of vector components
n_R   = 101 !< total number of radial points
n_Z   = 101 !< total number of vertical coordinate points
fields_filename = 'jorek_equilibrium'

!> Initialisation ---------------------------------------------------------------------------------------
!> initialise the MPI communicator
call init_mpi_threads(my_id,n_cpus,ierr)
!> read MHD fields
write(*,*) "Reading MHD data ..."
call sim%initialize(0,.true.,my_id,n_cpus)
field_reader = event(read_jorek_fields_interp_linear(basename=trim(fields_filename),i=-1))
call with(sim,field_reader)
write(*,*) "Reading MHD data: completed!"
!> initialise magnetic field array
write(*,*) "Initialise magnetic field datastructures ..."
!> allocate R,Z mesh and magnetic field arrays
allocate(R_mesh(n_R)); R_mesh = 0d0; allocate(Z_mesh(n_Z)); Z_mesh = 0d0;
allocate(magnetic_field(n_vec,n_Z,n_R)); magnetic_field = 0d0;
write(*,*) "Initialise magnetic field datastructures: completed!"

!> Compute mesh -----------------------------------------------------------------------------------------
write(*,*) "Generate computational mesh ..."
call generate_equidistant_RZ_mesh(sim%fields,n_R,n_Z,R_mesh,Z_mesh)
write(*,*) "Generate computational mesh: completed!"

!> Write input files ------------------------------------------------------------------------------------
#ifdef USE_HDF5
  write(*,*) "Write magnetic field in HDF5 file ..."
  write(*,*) "Write magnetic field in HDF5 file: completed!"
  write(*,*) "Write probability distribution function in HDF5 file ..."
  write(*,*) "Write probability distribution function in HDF5 file: completed!"
#endif
!> Finalisation -----------------------------------------------------------------------------------------
if(allocated(R_mesh))         deallocate(R_mesh);
if(allocated(Z_mesh))         deallocate(Z_mesh);
if(allocated(magnetic_field)) deallocate(magnetic_field);
call finalize_mpi_threads(ierr)

contains

!> Tools ------------------------------------------------------------------------------------------------
!> generate equidistant mesh in R and Z coordinates
!> inputs:
!>   fields: (fields_base) jorek mhd fields datatype 
!>   n_R:    (integer) number of major radius nodes
!>   n_Z:    (integer) number of vertical coordinate nodes
!> outputs:
!>   R_mesh: (real8)(n_R) major radius mesh
!>   Z_mesh: (real8)(n_Z) vertical coordinate mesh
subroutine generate_equidistant_RZ_mesh(fields,n_R,n_Z,R_mesh,Z_mesh)
use mod_fields, only: fields_base
implicit none
!> inputs:
class(fields_base),intent(in) :: fields
integer,intent(in)            :: n_R,n_Z
!> outputs: 
real*8,dimension(n_R),intent(out) :: R_mesh
real*8,dimension(n_Z),intent(out) :: Z_mesh
!> variables:
integer :: ii
real*8  :: dR,dZ
!> define the domain bounding box
call domain_bounding_box(fields%node_list,fields%element_list,R_mesh(1),R_mesh(n_R),Z_mesh(1),Z_mesh(n_Z))
dR = (R_mesh(n_R)-R_mesh(1))/real(n_R-1,kind=8); dZ = (Z_mesh(n_Z)-Z_mesh(1))/real(n_Z-1,kind=8);
write(*,*) "initial bounding box R,Z: ",R_mesh(1),R_mesh(n_R),Z_mesh(1),Z_mesh(n_Z)
!> generate mesh
R_mesh = [(R_mesh(1)+real(ii,kind=8)*dR,ii=0,n_R-1)]
Z_mesh = [(Z_mesh(1)+real(ii,kind=8)*dZ,ii=0,n_Z-1)]
write(*,*) "final bounding box R,Z: ",R_mesh(1),R_mesh(n_R),Z_mesh(1),Z_mesh(n_Z)
end subroutine generate_equidistant_RZ_mesh
!> ------------------------------------------------------------------------------------------------------
end program generate_input_files_for_SOFT


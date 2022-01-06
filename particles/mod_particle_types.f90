!> Mod_particle_type contains the default particle type, type_particle
!> If you wish to add a new particle type, create one here, inheriting from
!> [[particle_base]], and update [[mod_particle_io]] (search for `particle_kinetic`
!> and add your particle at each spot)
module mod_particle_types
  implicit none
  private
  public :: particle_base_id,particle_fieldline_id,particle_gc_id
  public :: particle_gc_vpar_id,particle_gc_Qin_id,particle_kinetic_id
  public :: particle_kinetic_leapfrog_id,particle_kinetic_relativistic_id
  public :: particle_gc_relativistic_id
  public :: particle_base, particle_kinetic, particle_kinetic_leapfrog
  public :: particle_gc, particle_fieldline
  public :: particle_kinetic_relativistic, particle_gc_relativistic
  public :: particle_gc_vpar, particle_gc_Qin
  public :: particle_get_q
  public :: copy_particle
  public :: copy_particle_base
  public :: copy_particle_kinetic_leapfrog
  public :: codify_particle_type
  public :: find_active_particle_id
  !> publicity only for unit testing
#ifdef UNIT_TESTS
  public :: find_active_particle_id_seq
  public :: find_active_particle_id_openmp
#endif

  !> enumerator of the particle type 
  enum, bind(C)
  enumerator :: particle_base_id=0,particle_fieldline_id,particle_gc_id,&
                particle_gc_vpar_id,particle_gc_Qin_id,particle_kinetic_id,&
                particle_kinetic_leapfrog_id,particle_kinetic_relativistic_id,&
                particle_gc_relativistic_id
  endenum

  !> The base type for all other particles. Includes only the position and weight elements
  !> Integration in a 2D finite element method is included in the form of 2 coordinates
  !> and an element index.
  type, abstract :: particle_base
    real*8    :: x(3)             !< particle position in real space
    real*8    :: st(2)            !< particle position in the element
    real*8    :: weight = 1.0     !< weight (i.e. number of particles)
    integer*4 :: i_elm = 0        !< index in element_list. Negative indices indicate lost particles on the edge of - that element.
    integer*4 :: i_life = 0       !< particle lifetime index (i.e. is this still the same particle?)
    real*4    :: t_birth = 0.0    !< birth time of this particle
    !< zero means lost without location specification.
  contains
    procedure :: copy => copy_particle
    generic :: assignment(=) => copy
  end type particle_base

  !> A simple type just for fieldline tracing in two-step methods (Adams Bashforth) or for forward euler
  type, extends(particle_base) :: particle_fieldline
    real*8    :: B_hat_prev(3) = 0.d0 !< Field direction at previous timestep
    real*8    :: v = 0.d0 !< Parallel velocity along the fieldline
  end type particle_fieldline

  !> A simple guiding-center particle type.
  type, extends(particle_base) :: particle_gc
    real*8    :: E  = 0.d0 !< The particle energy [eV]
    real*8    :: mu = 0.d0 !< The magnetic moment [eV/T]. Sign determines sign of v_par
    integer*1 :: q  = 0_1  !< Charge [e]
  end type particle_gc

  !> A simple guiding-center particle type.
  type, extends(particle_base) :: particle_gc_vpar
    real*8    :: vpar = 0.d0 !< Guiding centre parallel velocity [m/s]
    real*8    :: mu   = 0.d0 !< The magnetic moment [eV/T] 
    integer*1 :: q    = 0_1  !< Charge [e]
  end type particle_gc_vpar

  !> A simple guiding-center particle type.
  type, extends(particle_gc_vpar) :: particle_gc_Qin
    real*8    :: x_m(3)        !< position (previous step)
    real*8    :: vpar_m        !< parallel velocity (previous step)
    real*8    :: Astar_m(3)    !< A* (previous step)
    real*8    :: Astar_k(3)    !< A*  (current step)
    real*8    :: dAstar_k(3,3) !< dA* (current step)
    real*8    :: Bn_k          !< B   (amplitude, current step)
    real*8    :: dBn_k(3)      !< dB  (derivatives of Bn, current step)
    real*8    :: Bnorm_k(3)    !< normalised B (current step)
    real*8    :: E_k(3)        !< electric field (current step)
  end type particle_gc_Qin

  !> For most kinetic methods the velocity is required at time \(t\)
  type, extends(particle_base) :: particle_kinetic
    real*8, dimension(3) :: v = 0.d0 !< Velocity [m/s]
    integer*1            :: q = 0_1 !< charge [e]
  end type particle_kinetic

  !> Leapfrog methods define the particle velocity at time \(t^{n-1/2}\)
  !> and are therefore incompatible with normal kinetic methods (but a conversion
  !> function should not be too difficult)
  type, extends(particle_base) :: particle_kinetic_leapfrog
    real*8, dimension(3) :: v = 0.d0 !< Velocity [m/s] at t=t^(n-1/2) (where the position is known at t^n)
    integer*1            :: q = 0_1 !< charge [e]
  end type particle_kinetic_leapfrog

  !> This particle type is used for computing the full orbit trajectory
  !> of a relativistic particle. Particle position and momentum are given at time \(t\)
 type, extends(particle_base) :: particle_kinetic_relativistic
    real(kind=8),dimension(3) :: p !< Momentum in Cartesian coordinates (p_x,p_y,p_z) in [AMU*m/s]
    integer(kind=1)           :: q !< charge [e]
 end type particle_kinetic_relativistic

 !> This particle type is used for computing the guiding center trajectory of a
 !> relativistic particle. GC position and momentum are given at the time (\t\)
 type, extends(particle_base) :: particle_gc_relativistic
    real(kind=8), dimension(2) :: p  !< 1: parallel momentum [AMU m/s], 2: magnetic moment [(AMU*m**2)/(T*s**2)]
    integer(kind=1) :: q !< charge [e]
 end type particle_gc_relativistic

!> interfaces ------------------------------------------------------------------------------
interface codify_particle_type
  module procedure codify_single_particle_type
  module procedure codify_particle_list_alloc_type
end interface codify_particle_type

interface find_active_particle_id
  module procedure find_active_particle_id_base
  module procedure find_active_particle_id_type
end interface find_active_particle_id

contains
  !> Convenience function to obtain q if it exists, or 0 otherwise
  !> Here also because of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=82064
  !> which means that we cannot use the same derived type in too many modules which will be
  !> imported in the main program (roughly)
  pure function particle_get_q(in) result(q)
    class(particle_base), intent(in) :: in
    integer*1 :: q
    select type (p => in)
    type is (particle_kinetic)
      q = p%q
    type is (particle_kinetic_leapfrog)
      q = p%q
    type is (particle_gc)
      q = p%q
    type is (particle_gc_vpar)
      q = p%q
    type is (particle_gc_Qin)
      q = p%q
    type is (particle_kinetic_relativistic)
      q = p%q
    type is (particle_gc_relativistic)
      q=p%q
    class default
      q = 0
    end select
  end function particle_get_q
  !> Copy the base variables from one particle of class(particle_base) to another
  
  pure subroutine copy_particle_base(in, out)
    class(particle_base), intent(in)    :: in
    class(particle_base), intent(inout) :: out
    out%x      = in%x
    out%st     = in%st
    out%weight = in%weight
    out%i_elm  = in%i_elm
    out%i_life = in%i_life
    out%t_birth= in%t_birth
  end subroutine copy_particle_base

  !> Copy one particle of a type kinetic_leapfrog to another
  pure subroutine copy_particle_kinetic_leapfrog(in, out)
    type(particle_kinetic_leapfrog), intent(in)    :: in
    type(particle_kinetic_leapfrog), intent(inout) :: out
    out%x       = in%x
    out%st      = in%st
    out%weight  = in%weight
    out%i_elm   = in%i_elm
    out%i_life  = in%i_life
    out%t_birth = in%t_birth
    out%v       = in%v
    out%q       = in%q
  end subroutine copy_particle_kinetic_leapfrog

  !> Copy a descendant of particle_base into another descendant of particle_base
  !> as requested by the types of the input and output parameters.
  !> We do not transform any of the quantities between eachother, just copy them
  !> if the particles are of the same type.
  !> To do a proper transform we typically need additional parameters, like the
  !> mass or magnetic field. See [[mod_boris.f90]] for some examples.
  !> See https://stackoverflow.com/a/19082934
  subroutine copy_particle(particle_out, particle_in)
    class(particle_base), intent(out) :: particle_out !< Particle to copy attributes into
    class(particle_base), intent(in)  :: particle_in  !< Particle to copy attributes from

    particle_out%x        = particle_in%x
    particle_out%st       = particle_in%st
    particle_out%weight   = particle_in%weight
    particle_out%i_elm    = particle_in%i_elm
    particle_out%i_life   = particle_in%i_life  !< to be checked
    particle_out%t_birth  = particle_in%t_birth !< to be checked

    select type (p_out => particle_out)
    type is (particle_fieldline)
      select type (p_in => particle_in)
      type is (particle_fieldline) ! this is a straight copy, simple
        p_out%B_hat_prev = p_in%B_hat_prev
        p_out%v = p_in%v
      class default
        ! Maybe we should warn here instead of just putting zeros,
        ! or put nothing so the compiler can catch the uninitialized value
        p_out%B_hat_prev = [0.d0, 0.d0, 0.d0]
        p_out%v = 0.d0
      end select
    type is (particle_gc)
      select type (p_in => particle_in)
      type is (particle_gc)
        p_out%E  = p_in%E
        p_out%mu = p_in%mu
        p_out%q  = p_in%q
      class default
        p_out%E  = 0.d0
        p_out%mu = 0.d0
        p_out%q  = 0
      end select
    type is (particle_gc_vpar)
      select type (p_in => particle_in)
      type is (particle_gc_vpar)
        p_out%vpar = p_in%vpar
        p_out%mu   = p_in%mu
        p_out%q    = p_in%q
      class default
        p_out%vpar = 0.d0
        p_out%mu   = 0.d0
        p_out%q    = 0
      end select
    type is (particle_gc_Qin)
    select type (p_in => particle_in)
    type is (particle_gc_Qin)
      p_out%vpar     = p_in%vpar
      p_out%mu       = p_in%mu
      p_out%q        = p_in%q
      p_out%x_m      = p_in%x_m
      p_out%vpar_m   = p_in%vpar_m
      p_out%Astar_m  = p_in%Astar_m
      p_out%Astar_k  = p_in%Astar_k
      p_out%dAstar_k = p_in%dAstar_k
      p_out%Bn_k     = p_in%Bn_k
      p_out%dBn_k    = p_in%dBn_k
      p_out%Bnorm_k  = p_in%Bnorm_k
      p_out%E_k      = p_in%E_k
    class default
      p_out%vpar     = 0.d0
      p_out%mu       = 0.d0
      p_out%q        = 0
      p_out%x_m      = 0.d0
      p_out%vpar_m   = 0.d0
      p_out%Astar_m  = 0.d0
      p_out%Astar_k  = 0.d0
      p_out%dAstar_k = 0.d0
      p_out%Bn_k     = 0.d0
      p_out%dBn_k    = 0.d0
      p_out%Bnorm_k  = 0.d0
      p_out%E_k      = 0.d0
    end select
    type is (particle_kinetic)
      select type (p_in => particle_in)
      type is (particle_kinetic)
        p_out%v  = p_in%v
        p_out%q  = p_in%q
      class default
        p_out%v  = [0.d0, 0.d0, 0.d0]
        p_out%q  = 0
      end select
    type is (particle_kinetic_leapfrog)
      select type (p_in => particle_in)
      type is (particle_kinetic_leapfrog)
        p_out%v  = p_in%v
        p_out%q  = p_in%q
      class default
        ! the transformation from kinetic to kinetic_leapfrog could be done with a small error here
        p_out%v  = [0.d0, 0.d0, 0.d0]
        p_out%q  = 0
      end select
    type is (particle_kinetic_relativistic)
      select type (p_in => particle_in)
      type is (particle_kinetic_relativistic)
        p_out%p  = p_in%p
        p_out%q  = p_in%q
      class default
        p_out%p = [0.d0,0.d0,0.d0]
        p_out%q = 0
      end select     
     type is (particle_gc_relativistic)
       select type (p_in => particle_in)
       type is (particle_gc_relativistic)
         p_out%p = p_in%p
         p_out%q = p_in%q
     class default
        p_out%p  = [0.d0, 0.d0]
        p_out%q  = 0
      end select
    end select
  end subroutine copy_particle

  !> codify_single_particle_type returns a codified 
  !> particle type for a particle
  !> condification:
  !>   0 -> default
  !>   1 -> particle_fieldline
  !>   2 -> particle_gc
  !>   3 -> particle_gc_vpar
  !>   4 -> particle_gc_Qin
  !>   5 -> particle_kinetic
  !>   6 -> particle_kinetic_leapfrog
  !>   7 -> particle_kinetic_relativistic
  !>   8 -> particle_gc_relativistic
  function codify_single_particle_type(particle) result(p_type)
    implicit none
    class(particle_base),intent(in) :: particle
    integer :: p_type
    p_type = 0
    select type (p=>particle)
      type is (particle_fieldline)
      p_type = particle_fieldline_id
      type is (particle_gc)
      p_type = particle_gc_id
      type is (particle_gc_vpar)
      p_type = particle_gc_vpar_id
      type is (particle_gc_Qin)
      p_type = particle_gc_Qin_id
      type is (particle_kinetic)
      p_type = particle_kinetic_id
      type is (particle_kinetic_leapfrog)
      p_type = particle_kinetic_leapfrog_id
      type is (particle_kinetic_relativistic)
      p_type = particle_kinetic_relativistic_id
      type is (particle_gc_relativistic)
      p_type = particle_gc_relativistic_id
    end select
  end function codify_single_particle_type
 
  !> codify_particle_list_alloc_type returns a 
  !> codified particle type for a particle list
  !> condification:
  !>   1 -> particle_fieldline
  !>   2 -> particle_gc
  !>   3 -> particle_gc_vpar
  !>   4 -> particle_gc_Qin
  !>   5 -> particle_kinetic
  !>   6 -> particle_kinetic_leapfrog
  !>   7 -> particle_kinetic_relativistic
  !>   8 -> particle_gc_relativistic
  function codify_particle_list_alloc_type(particle_list) result(p_type)
    implicit none
    class(particle_base),dimension(:),allocatable,intent(in) :: particle_list
    integer :: p_type
    p_type = 0
    select type (p=>particle_list)
      type is (particle_fieldline)
      p_type = particle_fieldline_id
      type is (particle_gc)
      p_type = particle_gc_id
      type is (particle_gc_vpar)
      p_type = particle_gc_vpar_id
      type is (particle_gc_Qin)
      p_type = particle_gc_Qin_id
      type is (particle_kinetic)
      p_type = particle_kinetic_id
      type is (particle_kinetic_leapfrog)
      p_type = particle_kinetic_leapfrog_id
      type is (particle_kinetic_relativistic)
      p_type = particle_kinetic_relativistic_id
      type is (particle_gc_relativistic)
      p_type = particle_gc_relativistic_id
    end select
  end function codify_particle_list_alloc_type

  !> find_active_particle_id_type returns the number 
  !> and index of acitve particles withing a list for 
  !> are specific type of particle (encoded).
  !> Active particles are particles having i_elm>0.
  !> If particle code - particle type do not match
  !> returns 0 
  !> inputs:
  !>   particle_code: (integer) encoded particle type
  !>   n_particles:   (integer) number of particles
  !>   particle_list: (particle_base)(n_particles) particle list
  !> outputs:
  !>   n_active_particles: (integer) number of active particles
  !>   active_particle_id: (integer)(n_particles) active particle indices
  subroutine find_active_particle_id_type(particle_code,n_particles,&
  particle_list,n_active_particles,active_particle_id)
    implicit none
    !> inputs
    integer,intent(in) :: particle_code,n_particles
    class(particle_base),dimension(:),allocatable,intent(in) :: particle_list
    !> outputs
    integer,intent(out) :: n_active_particles
    integer,dimension(n_particles),intent(out) :: active_particle_id
    !> use select type as if condition, ok it is ugly ...
    n_active_particles = 0; active_particle_id = 0;
    select type (p_list=>particle_list)
      type is (particle_fieldline)
        if(particle_code.eq.particle_fieldline_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_gc)
        if(particle_code.eq.particle_gc_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_gc_vpar)
        if(particle_code.eq.particle_gc_vpar_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_gc_Qin)
        if(particle_code.eq.particle_gc_Qin_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_kinetic)
        if(particle_code.eq.particle_kinetic_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_kinetic_leapfrog)
        if(particle_code.eq.particle_kinetic_leapfrog_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_kinetic_relativistic)
        if(particle_code.eq.particle_kinetic_relativistic_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
      type is (particle_gc_relativistic)
        if(particle_code.eq.particle_gc_relativistic_id) &
          call find_active_particle_id_base(n_particles,particle_list,&
               n_active_particles,active_particle_id)
    end select
  end subroutine find_active_particle_id_type

  !> find_active_particle_id_base returns the number 
  !> and index of acitve particles withing a list
  !> active particles are particles having i_elm>0
  !> interface for the sequential and openmp versions
  !> inputs:
  !>   n_particles:   (integer) number of particles
  !>   particle_list: (particle_base)(n_particles) particle list
  !> outputs:
  !>   n_active_particles: (integer) number of active particles
  !>   active_particle_id: (integer)(n_particles) active particle indices
  subroutine find_active_particle_id_base(n_particles,particle_list,&
  n_active_particles,active_particle_id)
    implicit none
    !> inputs
    integer,intent(in) :: n_particles
    class(particle_base),dimension(:),allocatable,intent(in) :: particle_list
    !> outputs
    integer,intent(out) :: n_active_particles
    integer,dimension(n_particles),intent(out) :: active_particle_id
#ifdef _OPENMP
    call find_active_particle_id_openmp(n_particles,particle_list,&
    n_active_particles,active_particle_id)
#else
    call find_active_particle_id_seq(n_particles,particle_list,&
    n_active_particles,active_particle_id)
#endif
  end subroutine find_active_particle_id_base

  !> find_active_particle_id_seq returns the number 
  !> and index of acitve particles withing a list
  !> active particles are particles having i_elm>0
  !> sequential version
  !> inputs:
  !>   n_particles:   (integer) number of particles
  !>   particle_list: (particle_base)(n_particles) particle list
  !> outputs:
  !>   n_active_particles: (integer) number of active particles
  !>   active_particle_id: (integer)(n_particles) active particle indices
  subroutine find_active_particle_id_seq(n_particles,particle_list,&
  n_active_particles,active_particle_id)
    implicit none
    !> inputs
    integer,intent(in) :: n_particles
    class(particle_base),dimension(:),allocatable,intent(in) :: particle_list
    !> outputs
    integer,intent(out) :: n_active_particles
    integer,dimension(n_particles),intent(out) :: active_particle_id
    !> variables
    integer :: ii
    n_active_particles = 0; active_particle_id = 0;
    do ii=1,n_particles
      if(particle_list(ii)%i_elm.lt.1) cycle !< skip invalid particle
       n_active_particles = n_active_particles + 1
       active_particle_id(n_active_particles) = ii
    enddo
  end subroutine find_active_particle_id_seq
  !> find_active_particle_id_openmp returns the number 
  !> and index of acitve particles withing a list
  !> active particles are particles having i_elm>0
  !> openmp enabled version.
  !> Note: the proposed openmp version is suboptimal and
  !>       hence, it must be improved in future
  !> inputs:
  !>   n_particles:   (integer) number of particles
  !>   particle_list: (particle_base)(n_particles) particle list
  !> outputs:
  !>   n_active_particles: (integer) number of active particles
  !>   active_particle_id: (integer)(n_particles) active particle indices
  subroutine find_active_particle_id_openmp(n_particles,particle_list,&
  n_active_particles,active_particle_id)
    !$ use omp_lib
    implicit none
    !> inputs
    integer,intent(in) :: n_particles
    class(particle_base),dimension(:),allocatable,intent(in) :: particle_list
    !> outputs
    integer,intent(out) :: n_active_particles
    integer,dimension(n_particles),intent(out) :: active_particle_id
    !> variables
    integer :: ii,jj,start_id
    integer :: thread_id,max_num_thread,thread_num,n_particles_thread,thread_num_out
    integer,dimension(:),allocatable :: n_active_particle_thread
    integer,dimension(:,:),allocatable :: particle_id_thread
    max_num_thread = 1
    !> allocate thread private arrays
    !$ max_num_thread = omp_get_max_threads()
    allocate(n_active_particle_thread(max_num_thread)); 
   !> assume heuristically that the actual number of particles per thread is not
   !> larger than 25% of the number of particles per thread computed using max_num_thread
    allocate(particle_id_thread(floor(1.25d0*real(n_particles/&
    max_num_thread,kind=8))+1,max_num_thread));
    n_active_particles = 0; active_particle_id = 0;
    n_active_particle_thread = 0

    !> find active particles and store their index for each thread
    !$omp parallel default(private) shared(n_particles,particle_list,&
    !$omp n_active_particle_thread,particle_id_thread,thread_num_out)
    thread_id = 1; thread_num = 1;
    !$ thread_id = omp_get_thread_num() + 1
    !$ thread_num = omp_get_num_threads()
    n_particles_thread = n_particles/thread_num
    start_id = n_particles_thread*(thread_id-1)
    if(thread_id.eq.thread_num) n_particles_thread = n_particles - n_particles_thread*(thread_num-1)
    !$omp master
    thread_num_out = thread_num !< use additional variable because lastprivate not supported
    !$omp end master
    do ii=start_id+1,start_id+n_particles_thread
      if(particle_list(ii)%i_elm.lt.1) cycle !< skip invalid particle
      n_active_particle_thread(thread_id) = n_active_particle_thread(thread_id) + 1
      particle_id_thread(n_active_particle_thread(thread_id),thread_id) = ii
    enddo
    !$omp end parallel
    !> assemble the particle id array
    do ii=1,thread_num_out
      active_particle_id(n_active_particles+1:n_active_particles+n_active_particle_thread(ii)) = &
      particle_id_thread(1:n_active_particle_thread(ii),ii)
      n_active_particles = n_active_particles + n_active_particle_thread(ii)
    enddo

    !> cleanup
    deallocate(n_active_particle_thread)
    deallocate(particle_id_thread)
  end subroutine find_active_particle_id_openmp

end module mod_particle_types

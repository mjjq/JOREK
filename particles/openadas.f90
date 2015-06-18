!------------------------------------------------------------------------------------------------------
! module takes the OPEN-ADAS data to calculate the steady state (or time evolution) charge distribution
! and average charge state as a function of temperature
! to be done : radiation losses
!
!
! Guido Huijsmans, ITER, 20/01/2015
!------------------------------------------------------------------------------------------------------

module openadas

integer             :: n_z, n_states, n_d_ion, n_T_ion, n_d_rec, n_T_rec, n_d_rb, n_T_rb, n_d_lt, n_T_lt
real*8, allocatable :: density_ion(:),    temperature_ion(:),    ionisation(:,:,:)
real*8, allocatable :: density_rec(:),    temperature_rec(:),    recombination(:,:,:)
real*8, allocatable :: density_rb(:),     temperature_rb(:),     rad_rb(:,:,:)
real*8, allocatable :: density_lt(:),     temperature_lt(:),     rad_lt(:,:,:)

integer :: n_d_cor, n_T_cor
real*8, allocatable :: density_cor(:), temperature_cor(:), Z_cor(:,:), rad_cor(:,:)

contains

subroutine read_adas

implicit none

integer :: i, izmin, izmax

write(*,'(A)') '*********************************'
write(*,'(A)') '* Importing OpenAdas data       *'
write(*,'(A)') '*********************************'

open(10,file='ionisation.txt')

read(10,*)  n_z, n_d_ion, n_T_ion, izmin, izmax

allocate(density_ion(n_d_ion), temperature_ion(n_T_ion), ionisation(n_Z+1,n_d_ion,n_T_ion))

read(10,*)

read(10,*) density_ion(:)
read(10,*) temperature_ion(:)

ionisation = -30.d0
do i = 1, n_z
  read(10,*)
  read(10,*) ionisation(i,:,:)
enddo

close(10)

open(10,file='recombination.txt')

read(10,*)  n_z, n_d_rec, n_T_rec, izmin, izmax

allocate(density_rec(n_d_rec), temperature_rec(n_T_rec), recombination(n_Z+1,n_d_rec,n_T_rec))

read(10,*)
read(10,*) density_rec(:)
read(10,*) temperature_rec(:)

recombination = -30.d0
do i = 1, n_z
  read(10,*)
  read(10,*) recombination(i+1,:,:)
enddo

open(10,file='rad_rb.txt')

read(10,*)  n_z, n_d_rb, n_T_rb, izmin, izmax

allocate(density_rb(n_d_rb), temperature_rb(n_T_rb), rad_rb(n_Z+1,n_d_rb,n_T_rb))

read(10,*)

read(10,*) density_rb(:)
read(10,*) temperature_rb(:)

rad_rb = -30.d0
do i = 1, n_z
  read(10,*)
  read(10,*) rad_rb(i+1,:,:)
enddo

open(10,file='rad_lt.txt')

read(10,*)  n_z, n_d_lt, n_T_lt, izmin, izmax

allocate(density_lt(n_d_lt), temperature_lt(n_T_lt), rad_lt(n_Z+1,n_d_lt,n_T_lt))

read(10,*)
read(10,*) density_lt(:)
read(10,*) temperature_lt(:)

rad_lt = -30.d0
do i = 1, n_z
  read(10,*)
  read(10,*) rad_lt(i+1,:,:)
enddo

n_states = n_Z + 1

!density_ion     = 10**6 * 10**density_ion  ! [m^-3]
!density_rec     = 10**6 * 10**density_rec  ! [m^-3]
!density_rad     = 10**6 * 10**density_rad  ! [m^-3]
!temperature_ion = 10**temperature_ion      ! [eV]
!temperature_rec = 10**temperature_rec      ! [eV]
!temperature_rad = 10**temperature_rad      ! [eV]
!ionisation      = 10**ionisation           ! [s^-1]
!recombination   = 10**recombination        ! [s^-1]
!rad_db          = 10**rad_db               ! [s^-1]
!rad_lt          = 10**rad_lt               ! [s^-1]

close(10)

write(*,*) 'Done reading adas data'

end subroutine


subroutine ion_rec_rad(n_points,density,temperature,Z,ion_rate,rec_rate,rad_rate_rb,rad_rate_lt)
! linear interpolation of density and temperature (assumes equidistant coordinates)

implicit none

integer :: n_points
real*8  :: density(n_points)         ! log10 of density
real*8  :: temperature(n_points)     ! log10 of temperature
integer :: Z(n_points)               ! charge state
real*8  :: ion_rate(n_points)        ! ionisation rate [1/s]
real*8  :: rec_rate(n_points)        ! recombination rate [1/s]
real*8  :: rad_rate_rb(n_points)     ! radiated power recombombination and Brehmstrahlung
real*8  :: rad_rate_lt(n_points)     ! radiated power line emmission

real*8  :: ion_low, ion_high, rec_low, rec_high, rad_rb_low, rad_rb_high,  rad_lt_low, rad_lt_high
integer :: index_d_ion(n_points),    index_T_ion(n_points),    index_d_rec(n_points),    index_T_rec(n_points)
integer :: index_d_rb(n_points),     index_T_rb(n_points),     index_d_lt(n_points),     index_T_lt(n_points), index_Z(n_points)
integer :: i


index_d_ion = floor((density     - density_ion(1))     / (density_ion(n_d_ion)       - density_ion(1))     * n_d_ion)  + 1
index_T_ion = floor((temperature - temperature_ion(1)) / (temperature_ion(n_T_ion)   - temperature_ion(1)) * n_T_ion)  + 1
index_d_rec = floor((density     - density_rec(1))     / (density_rec(n_d_rec)       - density_rec(1))     * n_d_rec)  + 1
index_T_rec = floor((temperature - temperature_rec(1)) / (temperature_rec(n_T_rec)   - temperature_rec(1)) * n_T_rec)  + 1
index_d_rb  = floor((density     - density_rb(1))      / (density_rb(n_d_rb)         - density_rb(1))      * n_d_rb)  + 1
index_T_rb  = floor((temperature - temperature_rb(1))  / (temperature_rb(n_T_rb)     - temperature_rb(1))  * n_T_rb)  + 1
index_d_lt  = floor((density     - density_lt(1))      / (density_lt(n_d_lt)         - density_lt(1))      * n_d_lt)  + 1
index_T_lt  = floor((temperature - temperature_lt(1))  / (temperature_lt(n_T_lt)     - temperature_lt(1))  * n_T_lt)  + 1
index_Z     = Z + 1

do i=1, n_points

  ion_low = ionisation(index_Z(i),index_d_ion(i),index_T_ion(i)) &
          + (ionisation(index_Z(i),index_d_ion(i)+1,index_T_ion(i)) -  ionisation(index_Z(i),index_d_ion(i),index_T_ion(i))) &
          * (density(i) - density_ion(index_d_ion(i))) / (density_ion(index_d_ion(i)+1) - density_ion(index_d_ion(i)))

  ion_high = ionisation(index_Z(i),index_d_ion(i),index_T_ion(i)+1) &
          + (ionisation(index_Z(i),index_d_ion(i)+1,index_T_ion(i)+1) -  ionisation(index_Z(i),index_d_ion(i),index_T_ion(i)+1)) &
          * (density(i) - density_ion(index_d_ion(i))) / (density_ion(index_d_ion(i)+1) - density_ion(index_d_ion(i)))

  ion_rate(i) = ion_low + (temperature(i) - temperature_ion(index_T_ion(i))) / (temperature_ion(index_T_ion(i)+1) - temperature_ion(index_T_ion(i))) &
                        * (ion_high - ion_low)

  rec_low = recombination(index_Z(i),index_d_rec(i),index_T_rec(i)) &
          + (recombination(index_Z(i),index_d_rec(i)+1,index_T_rec(i)) -  recombination(index_Z(i),index_d_rec(i),index_T_rec(i))) &
          * (density(i) - density_rec(index_d_rec(i))) / (density_rec(index_d_rec(i)+1) - density_rec(index_d_rec(i)))

  rec_high = recombination(index_Z(i),index_d_rec(i),index_T_rec(i)+1) &
          + (recombination(index_Z(i),index_d_rec(i)+1,index_T_rec(i)+1) -  recombination(index_Z(i),index_d_rec(i),index_T_rec(i)+1)) &
          * (density(i) - density_rec(index_d_rec(i))) / (density_rec(index_d_rec(i)+1) - density_rec(index_d_rec(i)))

  rec_rate(i) = rec_low + (temperature(i) - temperature_rec(index_T_rec(i))) / (temperature_rec(index_T_rec(i)+1) - temperature_rec(index_T_rec(i))) &
                        * (rec_high - rec_low)

  rad_rb_low = rad_rb(index_Z(i),index_d_rb(i),index_T_rb(i)) &
             + (rad_rb(index_Z(i),index_d_rb(i)+1,index_T_rb(i)) -  rad_rb(index_Z(i),index_d_rb(i),index_T_rb(i))) &
             * (density(i) - density_rb(index_d_rb(i))) / (density_rb(index_d_rb(i)+1) - density_rb(index_d_rb(i)))

  rad_rb_high = rad_rb(index_Z(i),index_d_rb(i),index_T_rb(i)+1) &
              + (rad_rb(index_Z(i),index_d_rb(i)+1,index_T_rb(i)+1) -  rad_rb(index_Z(i),index_d_rb(i),index_T_rb(i)+1)) &
              * (density(i) - density_rb(index_d_rb(i))) / (density_rb(index_d_rb(i)+1) - density_rb(index_d_rb(i)))

  rad_rate_rb(i) = rad_rb_low + (temperature(i) - temperature_rb(index_T_rb(i))) / (temperature_rb(index_T_rb(i)+1) - temperature_rb(index_T_rb(i))) &
              * (rad_rb_high - rad_rb_low)

  rad_lt_low = rad_lt(index_Z(i),index_d_lt(i),index_T_lt(i)) &
             + (rad_lt(index_Z(i),index_d_lt(i)+1,index_T_lt(i)) -  rad_lt(index_Z(i),index_d_lt(i),index_T_lt(i))) &
             * (density(i) - density_lt(index_d_lt(i))) / (density_lt(index_d_lt(i)+1) - density_lt(index_d_lt(i)))

  rad_lt_high = rad_lt(index_Z(i),index_d_lt(i),index_T_lt(i)+1) &
              + (rad_lt(index_Z(i),index_d_lt(i)+1,index_T_lt(i)+1) -  rad_lt(index_Z(i),index_d_lt(i),index_T_lt(i)+1)) &
              * (density(i) - density_lt(index_d_lt(i))) / (density_lt(index_d_lt(i)+1) - density_lt(index_d_lt(i)))

  rad_rate_lt(i) = rad_lt_low + (temperature(i) - temperature_lt(index_T_lt(i))) / (temperature_lt(index_T_lt(i)+1) - temperature_lt(index_T_lt(i))) &
              * (rad_lt_high - rad_lt_low)

enddo

do i=1, n_points
  ion_rate(i)    = 10**ion_rate(i)
  rec_rate(i)    = 10**rec_rate(i)
  rad_rate_rb(i) = 10**rad_rate_rb(i)
  rad_rate_lt(i) = 10**rad_rate_lt(i)
enddo

return
end subroutine

subroutine interpolate_ion_rec_rad(density,temperature,ion_rate,rec_rate,rad_rate_rb,rad_rate_lt)
! linear interpolation of density and temperature

implicit none

real*8, allocatable  :: ion_rate(:), rec_rate(:), rad_rate_rb(:), rad_rate_lt(:)
real*8               :: density, temperature, ion_low, ion_high, rec_low, rec_high, rad_rb_low, rad_rb_high, rad_lt_low, rad_lt_high
integer              :: i, index_d_ion, index_T_ion, index_d_rec, index_T_rec, index_d_rb, index_T_rb, index_d_lt, index_T_lt

index_d_ion = minloc(abs(density-density_ion),1)
index_T_ion = minloc(abs(temperature-temperature_ion),1)
index_d_rec = minloc(abs(density-density_rec),1)
index_T_rec = minloc(abs(temperature-temperature_rec),1)
index_d_rb  = minloc(abs(density-density_rb),1)
index_T_rb  = minloc(abs(temperature-temperature_rb),1)
index_d_lt  = minloc(abs(density-density_lt),1)
index_T_lt  = minloc(abs(temperature-temperature_lt),1)

if (density_ion(index_d_ion)     - density     .ge. 0.d0) index_d_ion = max(1,index_d_ion - 1)
if (temperature_ion(index_T_ion) - temperature .ge. 0.d0) index_T_ion = max(1,index_T_ion - 1)
if (density_rec(index_d_rec)     - density     .ge. 0.d0) index_d_rec = max(1,index_d_rec - 1)
if (temperature_rec(index_T_rec) - temperature .ge. 0.d0) index_T_rec = max(1,index_T_rec - 1)
if (density_rb(index_d_rb)       - density     .ge. 0.d0) index_d_rb  = max(1,index_d_rb  - 1)
if (temperature_rb(index_T_rb)   - temperature .ge. 0.d0) index_T_rb  = max(1,index_T_rb  - 1)
if (density_lt(index_d_lt)       - density     .ge. 0.d0) index_d_lt  = max(1,index_d_lt  - 1)
if (temperature_lt(index_T_lt)   - temperature .ge. 0.d0) index_T_lt  = max(1,index_T_lt  - 1)

if (.not. allocated(ion_rate))    allocate(ion_rate(n_states))
if (.not. allocated(rec_rate))    allocate(rec_rate(n_states))
if (.not. allocated(rad_rate_rb)) allocate(rad_rate_rb(n_states))
if (.not. allocated(rad_rate_lt)) allocate(rad_rate_lt(n_states))

do i = 1, n_states

  ion_low = ionisation(i,index_d_ion,index_T_ion) &
          + (ionisation(i,index_d_ion+1,index_T_ion) -  ionisation(i,index_d_ion,index_T_ion)) &
          * (density - density_ion(index_d_ion)) / (density_ion(index_d_ion+1) - density_ion(index_d_ion))

  ion_high = ionisation(i,index_d_ion,index_T_ion+1) &
          + (ionisation(i,index_d_ion+1,index_T_ion+1) -  ionisation(i,index_d_ion,index_T_ion+1)) &
          * (density - density_ion(index_d_ion)) / (density_ion(index_d_ion+1) - density_ion(index_d_ion))

  ion_rate(i) = ion_low + (temperature - temperature_ion(index_T_ion)) / (temperature_ion(index_T_ion+1) -temperature_ion(index_T_ion)) &
                        * (ion_high - ion_low)

  rec_low = recombination(i,index_d_rec,index_T_rec) &
          + (recombination(i,index_d_rec+1,index_T_rec) -  recombination(i,index_d_rec,index_T_rec)) &
          * (density - density_rec(index_d_rec)) / (density_rec(index_d_rec+1) - density_rec(index_d_rec))

  rec_high = recombination(i,index_d_rec,index_T_rec+1) &
          + (recombination(i,index_d_rec+1,index_T_rec+1) -  recombination(i,index_d_rec,index_T_rec+1)) &
          * (density - density_rec(index_d_rec)) / (density_rec(index_d_rec+1) - density_rec(index_d_rec))

  rec_rate(i) = rec_low + (temperature - temperature_rec(index_T_rec)) / (temperature_rec(index_T_rec+1) -temperature_rec(index_T_rec)) &
                        * (rec_high - rec_low)

  rad_rb_low = rad_rb(i,index_d_rb,index_T_rb) &
          + (rad_rb(i,index_d_rb+1,index_T_rb) -  rad_rb(i,index_d_rb,index_T_rb)) &
          * (density - density_rb(index_d_rb)) / (density_rb(index_d_rb+1) - density_rb(index_d_rb))

  rad_rb_high = rad_rb(i,index_d_rb,index_T_rb+1) &
          + (rad_rb(i,index_d_rb+1,index_T_rb+1) - rad_rb(i,index_d_rb,index_T_rb+1)) &
          * (density - density_rb(index_d_rb)) / (density_rb(index_d_rb+1) - density_rb(index_d_rb))

  rad_rate_rb(i) = rad_rb_low + (temperature - temperature_rb(index_T_rb)) / (temperature_rb(index_T_rb+1) -temperature_rb(index_T_rb)) &
                 * (rad_rb_high - rad_rb_low)

  rad_lt_low = rad_lt(i,index_d_lt,index_T_lt) &
          + (rad_lt(i,index_d_lt+1,index_T_lt) -  rad_lt(i,index_d_lt,index_T_lt)) &
          * (density - density_lt(index_d_lt)) / (density_lt(index_d_lt+1) - density_lt(index_d_lt))

  rad_lt_high = rad_lt(i,index_d_lt,index_T_lt+1) &
          + (rad_lt(i,index_d_lt+1,index_T_lt+1) - rad_lt(i,index_d_lt,index_T_lt+1)) &
          * (density - density_lt(index_d_lt)) / (density_lt(index_d_lt+1) - density_lt(index_d_lt))

  rad_rate_lt(i) = rad_lt_low + (temperature - temperature_lt(index_T_lt)) / (temperature_lt(index_T_lt+1) -temperature_lt(index_T_lt)) &
                 * (rad_lt_high - rad_lt_low)

enddo

ion_rate    = 10**ion_rate          ! [s^-1]
rec_rate    = 10**rec_rate          ! [s^-1]
rad_rate_rb = 10**rad_rate_rb       ! [s^-1]
rad_rate_lt = 10**rad_rate_lt       ! [s^-1]

end subroutine

subroutine coronal

implicit none

real*8, allocatable :: coronal_Z(:),corona_matrix(:,:), zni(:), b(:), Z(:), fractions(:)
real*8, allocatable :: A_l(:),A_d(:),A_u(:),T_bar(:),Z_bar(:), R_bar(:)
real*8              :: density, temperature, log_density, log_temperature

real*8, allocatable :: ion_rate(:), rec_rate(:), rad_rate_rb(:), rad_rate_lt(:), Z_states(:)
real*8              :: t_step, theta, zc(30)
integer             :: n_states, i, info, n_step, j, k, m
integer             :: i_d_ion, i_T_ion, i_d_rec, i_T_rec

n_states = n_Z + 1

write(*,'(A)')      '*********************************'
write(*,'(A,i3,A)') '* Coronal model : ',n_states,'           *'
write(*,'(A)')      '*********************************'

n_d_cor = 4
n_T_cor = 200

allocate(density_cor(n_d_cor), temperature_cor(n_T_cor), Z_cor(n_d_cor,n_T_cor), rad_cor(n_d_cor,n_T_cor))

allocate(corona_matrix(3,n_states),zni(n_states),b(n_states), Z(n_states), fractions(n_states))

allocate(A_l(n_states-1),A_d(n_states),A_u(n_states-1))

t_step = 1.d6
theta  = 1.

density_cor     = (/ 12., 13., 14., 15. /)                                                             ! [m^3]

do m=1, n_d_cor

  density = 10**density_cor(m)

  do k=1, n_T_cor

    temperature_cor(k) = alog10( 1.d0 + exp(alog(4.d4 - 1.d0 + 1.d0)*float(k-1)/(float(n_T_cor-1))) - 1.d0 )  ! in [K]

    call interpolate_ion_rec_rad(density_cor(m),temperature_cor(k),ion_rate,rec_rate,rad_rate_rb,rad_rate_lt)

    corona_matrix = 0.d0

    do i=2, n_Z
      Z(i) = float(i-1)
      corona_matrix(1,i) = + ion_rate(i-1)
      corona_matrix(2,i) = - ion_rate(i)   - rec_rate(i)
      corona_matrix(3,i) =                 + rec_rate(i+1)
    enddo
    Z(1)        = 0.
    Z(n_states) = float(n_Z)

    corona_matrix(1,1)     = 0.
    corona_matrix(2,1)     = - ion_rate(1)      ! ionisation losses from n=0 to 1
    corona_matrix(3,1)     = + rec_rate(2)      ! recombination from n=1 to 0

    corona_matrix(1,n_Z+1) = + ion_rate(n_Z)
    corona_matrix(2,n_Z+1) = - rec_rate(n_Z+1)
    corona_matrix(3,n_Z+1) = 0.

    corona_matrix = corona_matrix * density

    zni = 1.e9

    do i=2,n_Z
      b(i) =  zni(i-1)*corona_matrix(1,i) + zni(i)*corona_matrix(2,i) + zni(i+1)*corona_matrix(3,i)
    enddo
    b(1)        = zni(1)*corona_matrix(2,1)                 + zni(2)*corona_matrix(3,1)
    b(n_states) = zni(n_states-1)*corona_matrix(1,n_states) + zni(n_states)*corona_matrix(2,n_states)

    A_l  =    - theta*t_step*corona_matrix(1,2:n_states)
    A_d  = 1. - theta*t_step*corona_matrix(2,:)
    A_u  =    - theta*t_step*corona_matrix(3,1:n_states-1)

    b = b * t_step

    call dgtsv(n_states,1,A_l,A_d,A_u,b,n_states,info)

    zni = zni + b

    if (info .ne. 0) write(*,*) 'info : ',info

    fractions    = zni / sum(zni)

    Z_cor(m,k)   = dot_product(fractions,Z)

    rad_cor(m,k) = alog10( dot_product(fractions,rad_rate_rb+rad_rate_lt))

  enddo

enddo

endsubroutine coronal

subroutine interpolate_coronal(density,temperature,Z_coronal,radiation_coronal)
! linear interpolation of density and temperature (assumes equidistant coordinates)

implicit none

real*8  :: density            ! log10 of density       [cm^-3]
real*8  :: temperature        ! log10 of temperature   [K]
real*8  :: Z_coronal          ! charge state of coronal equilibrium
real*8  :: radiation_coronal  ! radiated power of coronal equilibrium [Wcm^3] (CHECK!)

real*8  :: zcor_low, zcor_high, rad_low, rad_high
integer :: index_d_cor,   index_T_cor

index_d_cor = floor((density     - density_cor(1))     / (density_cor(n_d_cor)       - density_cor(1))     * n_d_cor)  + 1
index_T_cor = floor((temperature - temperature_cor(1)) / (temperature_cor(n_T_cor)   - temperature_cor(1)) * n_T_cor)  + 1

index_d_cor = max(1,min(index_d_cor,n_d_cor-1))
index_T_cor = max(1,min(index_T_cor,n_T_cor-1))

Zcor_low =   Z_cor(index_d_cor,index_T_cor) &
          + (Z_cor(index_d_cor+1,index_T_cor) -  Z_cor(index_d_cor,index_T_cor)) &
          * (density - density_cor(index_d_cor)) / (density_cor(index_d_cor+1) - density_cor(index_d_cor))

Zcor_high = Z_cor(index_d_cor,index_T_cor+1) &
          + (Z_cor(index_d_cor+1,index_T_cor+1) -  Z_cor(index_d_cor,index_T_cor+1)) &
          * (density - density_cor(index_d_cor)) / (density_cor(index_d_cor+1) - density_cor(index_d_cor))

Z_coronal = Zcor_low + (temperature - temperature_cor(index_T_cor)) / (temperature_cor(index_T_cor+1) - temperature_cor(index_T_cor)) &
                    * (Zcor_high - Zcor_low)

rad_low =    rad_cor(index_d_cor,index_T_cor) &
          + (rad_cor(index_d_cor+1,index_T_cor) -  rad_cor(index_d_cor,index_T_cor)) &
          * (density - density_cor(index_d_cor)) / (density_cor(index_d_cor+1) - density_cor(index_d_cor))

rad_high =   rad_cor(index_d_cor,index_T_cor+1) &
          + (rad_cor(index_d_cor+1,index_T_cor+1) -  rad_cor(index_d_cor,index_T_cor+1)) &
          * (density - density_cor(index_d_cor)) / (density_cor(index_d_cor+1) - density_cor(index_d_cor))

radiation_coronal = rad_low + (temperature - temperature_cor(index_T_cor)) / (temperature_cor(index_T_cor+1) - temperature_cor(index_T_cor)) &
                        * (rad_high - rad_low)

radiation_coronal = 10**radiation_coronal

return
endsubroutine interpolate_coronal

end module openadas


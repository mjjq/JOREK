!> Given (psi,R,Z)-derivatives of a variable, project on the node values
subroutine project_var_on_node(node_list, i_node, i_var, F_values)

use data_structure
use mod_basisfunctions
use mod_parameters, only: n_degrees, n_order
implicit none

! --- Routine parameters
type (type_node_list),    intent(inout):: node_list                    !< node list, duh
integer,                  intent(in)   :: i_node                       !< index of node to project on
integer,                  intent(in)   :: i_var                        !< variable we are projecting one
real*8,                   intent(in)   :: F_values(0:4,0:4,0:4)        !< variable and its derivatives
             !< The 3 indices are for the thre variables (psi, R, Z)
             !< most variables only use psi and Z, but for the current in equilibrium.f90
             !< there is also a dependence on R. Hence we just do the generic case for
             !< all options. The indices go from 0:4, 0 means no derivative, 4 means 4th
             !< derivative. So F_values(0,0,0) is the value at the node, F(0,1,0) is
             !< dF_dR, and F(2,1,1) is dF_dpsi2_dR_dZ.

! --- Variable derivatives
real*8  :: F
real*8  :: dF_dpsi, dF_dR, dF_dZ                                                         ! 1st derivatives
real*8  :: dF_dpsi2, dF_dpsi_dR, dF_dpsi_dZ                                              ! 2nd derivatives
real*8  :: dF_dR2, dF_dR_dZ, dF_dZ2                                                      ! 2nd derivatives
real*8  :: dF_dpsi3, dF_dpsi2_dR, dF_dpsi2_dZ, dF_dpsi_dR2, dF_dpsi_dZ2, dF_dpsi_dR_dZ   ! 3rd derivatives
real*8  :: dF_dR3, dF_dR2_dZ, dF_dR_dZ2, dF_dZ3                                          ! 3rd derivatives
real*8  :: dF_dpsi4, dF_dpsi3_dR, dF_dpsi3_dZ, dF_dpsi2_dR2, dF_dpsi2_dZ2, dF_dpsi2_dR_dZ! 4th derivatives
real*8  :: dF_dpsi_dR3, dF_dpsi_dR2_dZ, dF_dpsi_dR_dZ2, dF_dpsi_dZ3                      ! 4th derivatives
real*8  :: dF_dR4, dF_dR3_dZ, dF_dR2_dZ2, dF_dR_dZ3, dF_dZ4                              ! 4th derivatives

! --- Output variables
real*8  :: var_out(n_degrees)

! --- Local element psi variables
real*8  :: psi
real*8  :: psi_s, psi_t           ! 1st derivatives
real*8  :: psi_ss, psi_tt, psi_st ! 2nd derivatives
real*8  :: psi_sst, psi_stt       ! 3rd derivatives
real*8  :: psi_sstt               ! 4th derivatives

! --- Local element spatial derivatives
real*8  :: R, Z
real*8  :: R_s, R_t, Z_s, Z_t                  ! 1st derivatives
real*8  :: R_ss, R_tt, R_st, Z_ss, Z_tt, Z_st  ! 2nd derivatives
real*8  :: R_sst, R_stt, Z_sst, Z_stt          ! 3rd derivatives
real*8  :: R_sstt, Z_sstt                      ! 4th derivatives


! --- Allocate user's variables
F = F_values(0,0,0)
dF_dpsi      = F_values(1,0,0); dF_dR          = F_values(0,1,0); dF_dZ          = F_values(0,0,1)
dF_dpsi2     = F_values(2,0,0); dF_dpsi_dR     = F_values(1,1,0); dF_dpsi_dZ     = F_values(1,0,1)
dF_dR2       = F_values(0,2,0); dF_dR_dZ       = F_values(0,1,1); dF_dZ2         = F_values(0,0,2)
dF_dpsi3     = F_values(3,0,0); dF_dpsi2_dR    = F_values(2,1,0); dF_dpsi2_dZ    = F_values(2,0,1)
dF_dpsi_dR2  = F_values(1,2,0); dF_dpsi_dZ2    = F_values(1,0,2); dF_dpsi_dR_dZ  = F_values(1,1,1)
dF_dR3       = F_values(0,3,0); dF_dR2_dZ      = F_values(0,2,1); dF_dR_dZ2      = F_values(0,1,2); dF_dZ3      = F_values(0,0,3)
dF_dpsi4     = F_values(4,0,0); dF_dpsi3_dR    = F_values(3,1,0); dF_dpsi3_dZ    = F_values(3,0,1)
dF_dpsi2_dR2 = F_values(2,2,0); dF_dpsi2_dZ2   = F_values(2,0,2); dF_dpsi2_dR_dZ = F_values(2,1,1)
dF_dpsi_dR3  = F_values(1,3,0); dF_dpsi_dR2_dZ = F_values(1,2,1); dF_dpsi_dR_dZ2 = F_values(1,1,2); dF_dpsi_dZ3 = F_values(1,0,3)
dF_dR4       = F_values(0,4,0); dF_dR3_dZ      = F_values(0,3,1); dF_dR2_dZ2     = F_values(0,2,2); dF_dR_dZ3   = F_values(0,1,3)
dF_dZ4       = F_values(0,0,4)

! --- Allocate psi variables (in doubt, please refer to definition of derivatives for degrees 7,8,9)
psi      = node_list%node(i_node)%values(1,1,var_psi)
psi_s    = node_list%node(i_node)%values(1,2,var_psi)
psi_t    = node_list%node(i_node)%values(1,3,var_psi)
psi_st   = node_list%node(i_node)%values(1,4,var_psi)
if (n_order .eq. 5) then
psi_ss   = node_list%node(i_node)%values(1,5,var_psi)
psi_tt   = node_list%node(i_node)%values(1,6,var_psi)
psi_sst  = node_list%node(i_node)%values(1,7,var_psi) - node_list%node(i_node)%values(1,3,var_psi)
psi_stt  = node_list%node(i_node)%values(1,8,var_psi) - node_list%node(i_node)%values(1,2,var_psi)
psi_sstt = node_list%node(i_node)%values(1,9,var_psi) - node_list%node(i_node)%values(1,5,var_psi) - node_list%node(i_node)%values(1,6,var_psi)
endif

! --- Allocate RZ variables
R      = node_list%node(i_node)%x(1,1,1)
R_s    = node_list%node(i_node)%x(1,2,1)
R_t    = node_list%node(i_node)%x(1,3,1)
R_st   = node_list%node(i_node)%x(1,4,1)
if (n_order .eq. 5) then
R_ss   = node_list%node(i_node)%x(1,5,1)
R_tt   = node_list%node(i_node)%x(1,6,1)
R_sst  = node_list%node(i_node)%x(1,7,1) - node_list%node(i_node)%x(1,3,1)
R_stt  = node_list%node(i_node)%x(1,8,1) - node_list%node(i_node)%x(1,2,1)
R_sstt = node_list%node(i_node)%x(1,9,1) - node_list%node(i_node)%x(1,5,1) - node_list%node(i_node)%x(1,6,1)
endif
Z      = node_list%node(i_node)%x(1,1,2)
Z_s    = node_list%node(i_node)%x(1,2,2)
Z_t    = node_list%node(i_node)%x(1,3,2)
Z_st   = node_list%node(i_node)%x(1,4,2)
if (n_order .eq. 5) then
Z_ss   = node_list%node(i_node)%x(1,5,2)
Z_tt   = node_list%node(i_node)%x(1,6,2)
Z_sst  = node_list%node(i_node)%x(1,7,2) - node_list%node(i_node)%x(1,3,2)
Z_stt  = node_list%node(i_node)%x(1,8,2) - node_list%node(i_node)%x(1,2,2)
Z_sstt = node_list%node(i_node)%x(1,9,2) - node_list%node(i_node)%x(1,5,2) - node_list%node(i_node)%x(1,6,2)
endif


! --- Initialise variable to zero
var_out = 0.d0


! --- 1st degree (value)
var_out(1) = F

! --- 2nd degree (s-derivative)
var_out(2) = (dF_dpsi * psi_s) + (dF_dR * R_s) + (dF_dZ * Z_s)

! --- 3rd degree (t-derivative)
var_out(3) = (dF_dpsi * psi_t) + (dF_dR * R_t) + (dF_dZ * Z_t)

! --- 4th degree (st-derivative)
var_out(4) =                                                                    &
      ( (dF_dpsi2 * psi_t) + (dF_dpsi_dR * R_t) + (dF_dpsi_dZ * Z_t) ) * psi_s  &
    +                             dF_dpsi                              * psi_st &
    + ( (dF_dpsi_dR * psi_t) + (dF_dR2 * R_t)   + (dF_dR_dZ * Z_t)   ) * R_s    &
    +                              dF_dR                               * R_st   &
    + ( (dF_dpsi_dZ * psi_t) + (dF_dR_dZ * R_t) + (dF_dZ2 * Z_t)     ) * Z_s    &
    +                              dF_dZ                               * Z_st

! --- Beyond this point, only quintic matters
if (n_order .eq. 5) then

! --- 5th degree (ss-derivative)
var_out(5) =                                                                    &
      ( (dF_dpsi2 * psi_s) + (dF_dpsi_dR * R_s) + (dF_dpsi_dZ * Z_s) ) * psi_s  &
    +                             dF_dpsi                              * psi_ss &
    + ( (dF_dpsi_dR * psi_s) + (dF_dR2 * R_s)   + (dF_dR_dZ * Z_s)   ) * R_s    &
    +                              dF_dR                               * R_ss   &
    + ( (dF_dpsi_dZ * psi_s) + (dF_dR_dZ * R_s) + (dF_dZ2 * Z_s)     ) * Z_s    &
    +                              dF_dZ                               * Z_ss

! --- 6th degree (tt-derivative)
var_out(6) =                                                                    &
      ( (dF_dpsi2 * psi_t) + (dF_dpsi_dR * R_t) + (dF_dpsi_dZ * Z_t) ) * psi_t  &
    +                             dF_dpsi                              * psi_tt &
    + ( (dF_dpsi_dR * psi_t) + (dF_dR2 * R_t)   + (dF_dR_dZ * Z_t)   ) * R_t    &
    +                              dF_dR                               * R_tt   &
    + ( (dF_dpsi_dZ * psi_t) + (dF_dR_dZ * R_t) + (dF_dZ2 * Z_t)     ) * Z_t    &
    +                              dF_dZ                               * Z_tt

! --- 7th degree (sst-derivative)
var_out(7) =                                                                                         &
    + ((dF_dpsi3 * psi_t) + (dF_dpsi2_dR * R_t) + (dF_dpsi2_dZ * Z_t))            * psi_s  * psi_s   &
    +                           dF_dpsi2                                    * 2.0 * psi_st * psi_s   &
    + ((dF_dpsi2_dR * psi_t) + (dF_dpsi_dR2 * R_t) + (dF_dpsi_dR_dZ * Z_t)) * 2.0 * R_s    * psi_s   &
    +                           dF_dpsi_dR                                  * 2.0 * R_st   * psi_s   &
    +                           dF_dpsi_dR                                  * 2.0 * R_s    * psi_st  &
    + ((dF_dpsi2_dZ * psi_t) + (dF_dpsi_dR_dZ * R_t) + (dF_dpsi_dZ2 * Z_t)) * 2.0 * Z_s    * psi_s   &
    +                           dF_dpsi_dZ                                  * 2.0 * Z_st   * psi_s   &
    +                           dF_dpsi_dZ                                  * 2.0 * Z_s    * psi_st  &
    + ((dF_dpsi2 * psi_t) + (dF_dpsi_dR * R_t) + (dF_dpsi_dZ * Z_t))                       * psi_ss  &
    +                           dF_dpsi                                                    * psi_sst &
    + ((dF_dpsi_dR2 * psi_t) + (dF_dR3 * R_t) + (dF_dR2_dZ * Z_t))                * R_s    * R_s     &
    +                           dF_dR2                                      * 2.0 * R_st   * R_s     &
    + ((dF_dpsi_dR_dZ * psi_t) + (dF_dR2_dZ * R_t) + (dF_dR_dZ2 * Z_t))     * 2.0 * Z_s    * R_s     &
    +                           dF_dR_dZ                                    * 2.0 * Z_st   * R_s     &
    +                           dF_dR_dZ                                    * 2.0 * Z_s    * R_st    &
    + ((dF_dpsi_dR * psi_t) + (dF_dR2 * R_t) + (dF_dR_dZ * Z_t))                           * R_ss    &
    +                           dF_dR                                                      * R_sst   &
    + ((dF_dpsi_dZ2 * psi_t) + (dF_dR_dZ2 * R_t) + (dF_dZ3 * Z_t))                * Z_s    * Z_s     &
    +                           dF_dZ2                                      * 2.0 * Z_st   * Z_s     &
    + ((dF_dpsi_dZ * psi_t) + (dF_dR_dZ * R_t) + (dF_dZ2 * Z_t))                           * Z_ss    &
    +                           dF_dZ                                                      * Z_sst   &
    + var_out(3)


! --- 8th degree (stt-derivative)
var_out(8) =                                                                                         &
    + ((dF_dpsi3 * psi_s) + (dF_dpsi2_dR * R_s) + (dF_dpsi2_dZ * Z_s))            * psi_t  * psi_t   &
    +                           dF_dpsi2                                    * 2.0 * psi_st * psi_t   &
    + ((dF_dpsi2_dR * psi_s) + (dF_dpsi_dR2 * R_s) + (dF_dpsi_dR_dZ * Z_s)) * 2.0 * R_t    * psi_t   &
    +                           dF_dpsi_dR                                  * 2.0 * R_st   * psi_t   &
    +                           dF_dpsi_dR                                  * 2.0 * R_t    * psi_st  &
    + ((dF_dpsi2_dZ * psi_s) + (dF_dpsi_dR_dZ * R_s) + (dF_dpsi_dZ2 * Z_s)) * 2.0 * Z_t    * psi_t   &
    +                           dF_dpsi_dZ                                  * 2.0 * Z_st   * psi_t   &
    +                           dF_dpsi_dZ                                  * 2.0 * Z_t    * psi_st  &
    + ((dF_dpsi2 * psi_s) + (dF_dpsi_dR * R_s) + (dF_dpsi_dZ * Z_s))                       * psi_tt  &
    +                           dF_dpsi                                                    * psi_stt &
    + ((dF_dpsi_dR2 * psi_s) + (dF_dR3 * R_s) + (dF_dR2_dZ * Z_s))                * R_t    * R_t     &
    +                           dF_dR2                                      * 2.0 * R_st   * R_t     &
    + ((dF_dpsi_dR_dZ * psi_s) + (dF_dR2_dZ * R_s) + (dF_dR_dZ2 * Z_s))     * 2.0 * Z_t    * R_t     &
    +                           dF_dR_dZ                                    * 2.0 * Z_st   * R_t     &
    +                           dF_dR_dZ                                    * 2.0 * Z_t    * R_st    &
    + ((dF_dpsi_dR * psi_s) + (dF_dR2 * R_s) + (dF_dR_dZ * Z_s))                           * R_tt    &
    +                           dF_dR                                                      * R_stt   &
    + ((dF_dpsi_dZ2 * psi_s) + (dF_dR_dZ2 * R_s) + (dF_dZ3 * Z_s))                * Z_t    * Z_t     &
    +                           dF_dZ2                                      * 2.0 * Z_st   * Z_t     &
    + ((dF_dpsi_dZ * psi_s) + (dF_dR_dZ * R_s) + (dF_dZ2 * Z_s))                           * Z_tt    &
    +                           dF_dZ                                                      * Z_stt   &
    + var_out(2)


! --- 9th degree (sstt-derivative) that's the fun one...
var_out(9) =                                                                                               &
    + ((dF_dpsi4 * psi_s) + (dF_dpsi3_dR * R_s) + (dF_dpsi3_dZ * Z_s))         * psi_s   * psi_t  * psi_t  &
    +                           dF_dpsi3                                       * psi_ss  * psi_t  * psi_t  &
    +                           dF_dpsi3                                 * 2.0 * psi_s   * psi_st * psi_t  &
    + ((dF_dpsi_dR3 * psi_s) + (dF_dR4 * R_s) + (dF_dR3_dZ * Z_s))             * R_s     * R_t    * R_t    &
    +                           dF_dR3                                         * R_ss    * R_t    * R_t    &
    +                           dF_dR3                                   * 2.0 * R_s     * R_st   * R_t    &
    + ((dF_dpsi_dZ3 * psi_s) + (dF_dR_dZ3 * R_s) + (dF_dZ4 * Z_s))             * Z_s     * Z_t    * Z_t    &
    +                           dF_dZ3                                         * Z_ss    * Z_t    * Z_t    &
    +                           dF_dZ3                                   * 2.0 * Z_s     * Z_st   * Z_t    &
    + ((dF_dpsi3_dR * psi_s) + (dF_dpsi2_dR2 * R_s) + (dF_dpsi2_dR_dZ * Z_s))  * (R_s    * psi_t  * psi_t  + 2.0 * psi_s  * R_t  * psi_t ) &
    +                           dF_dpsi2_dR                                    * (R_ss   * psi_t  * psi_t  + 2.0 * psi_ss * R_t  * psi_t ) &
    +                           dF_dpsi2_dR                                    * (R_s    * psi_st * psi_t  + 2.0 * psi_s  * R_st * psi_t ) &
    +                           dF_dpsi2_dR                                    * (R_s    * psi_t  * psi_st + 2.0 * psi_s  * R_t  * psi_st) &
    + ((dF_dpsi3_dZ * psi_s) + (dF_dpsi2_dR_dZ * R_s) + (dF_dpsi2_dZ2 * Z_s))  * (Z_s    * psi_t  * psi_t  + 2.0 * psi_s  * Z_t  * psi_t ) &
    +                           dF_dpsi2_dZ                                    * (Z_ss   * psi_t  * psi_t  + 2.0 * psi_ss * Z_t  * psi_t ) &
    +                           dF_dpsi2_dZ                                    * (Z_s    * psi_st * psi_t  + 2.0 * psi_s  * Z_st * psi_t ) &
    +                           dF_dpsi2_dZ                                    * (Z_s    * psi_t  * psi_st + 2.0 * psi_s  * Z_t  * psi_st) &
    + ((dF_dpsi2_dR2 * psi_s) + (dF_dpsi_dR3 * R_s) + (dF_dpsi_dR2_dZ * Z_s))  * (psi_s  * R_t    * R_t    + 2.0 * R_s    * R_t  * psi_t ) &
    +                           dF_dpsi_dR2                                    * (psi_ss * R_t    * R_t    + 2.0 * R_ss   * R_t  * psi_t ) &
    +                           dF_dpsi_dR2                                    * (psi_s  * R_st   * R_t    + 2.0 * R_s    * R_st * psi_t ) &
    +                           dF_dpsi_dR2                                    * (psi_s  * R_t    * R_st   + 2.0 * R_s    * R_t  * psi_st) &
    + ((dF_dpsi2_dZ2 * psi_s) + (dF_dpsi_dR_dZ2 * R_s) + (dF_dpsi_dZ3 * Z_s))  * (psi_s  * Z_t    * Z_t    + 2.0 * Z_s    * Z_t  * psi_t ) &
    +                           dF_dpsi_dZ2                                    * (psi_ss * Z_t    * Z_t    + 2.0 * Z_ss   * Z_t  * psi_t ) &
    +                           dF_dpsi_dZ2                                    * (psi_s  * Z_st   * Z_t    + 2.0 * Z_s    * Z_st * psi_t ) &
    +                           dF_dpsi_dZ2                                    * (psi_s  * Z_t    * Z_st   + 2.0 * Z_s    * Z_t  * psi_st) &
    + ((dF_dpsi_dR2_dZ * psi_s) + (dF_dR3_dZ * R_s) + (dF_dR2_dZ2 * Z_s))      * (Z_s    * R_t    * R_t    + 2.0 * R_s    * Z_t  * R_t   ) &
    +                           dF_dR2_dZ                                      * (Z_ss   * R_t    * R_t    + 2.0 * R_ss   * Z_t  * R_t   ) &
    +                           dF_dR2_dZ                                      * (Z_s    * R_st   * R_t    + 2.0 * R_s    * Z_st * R_t   ) &
    +                           dF_dR2_dZ                                      * (Z_s    * R_t    * R_st   + 2.0 * R_s    * Z_t  * R_st  ) &
    + ((dF_dpsi_dR_dZ2 * psi_s) + (dF_dR2_dZ2 * R_s) + (dF_dR_dZ3 * Z_s))      * (R_s    * Z_t    * Z_t    + 2.0 * Z_s    * Z_t  * R_t   ) &
    +                           dF_dR_dZ2                                      * (R_ss   * Z_t    * Z_t    + 2.0 * Z_ss   * Z_t  * R_t   ) &
    +                           dF_dR_dZ2                                      * (R_s    * Z_st   * Z_t    + 2.0 * Z_s    * Z_st * R_t   ) &
    +                           dF_dR_dZ2                                      * (R_s    * Z_t    * Z_st   + 2.0 * Z_s    * Z_t  * R_st  ) &
    + ((dF_dpsi2_dR_dZ * psi_s) + (dF_dpsi_dR2_dZ * R_s) + (dF_dpsi_dR_dZ2 * Z_s)) &
                                                              * 2.0 * (Z_s  * R_t  * psi_t  + R_s  * Z_t  * psi_t  + psi_s  * Z_t  * R_t ) &
    +                           dF_dpsi_dR_dZ                 * 2.0 * (Z_ss * R_t  * psi_t  + R_ss * Z_t  * psi_t  + psi_ss * Z_t  * R_t ) &
    +                           dF_dpsi_dR_dZ                 * 2.0 * (Z_s  * R_st * psi_t  + R_s  * Z_st * psi_t  + psi_s  * Z_st * R_t ) &
    +                           dF_dpsi_dR_dZ                 * 2.0 * (Z_s  * R_t  * psi_st + R_s  * Z_t  * psi_st + psi_s  * Z_t  * R_st) &
    + ((dF_dpsi3 * psi_s) + (dF_dpsi2_dR * R_s) + (dF_dpsi2_dZ * Z_s))         * (psi_s  * psi_tt  + 2.0 * psi_st  * psi_t ) &
    +                           dF_dpsi2                                       * (psi_ss * psi_tt  + 2.0 * psi_sst * psi_t ) &
    +                           dF_dpsi2                                       * (psi_s  * psi_stt + 2.0 * psi_st  * psi_st) &
    + ((dF_dpsi_dR2 * psi_s) + (dF_dR3 * R_s) + (dF_dR2_dZ * Z_s))             * (R_s    * R_tt    + 2.0 * R_st    * R_t   ) &
    +                           dF_dR2                                         * (R_ss   * R_tt    + 2.0 * R_sst   * R_t   ) &
    +                           dF_dR2                                         * (R_s    * R_stt   + 2.0 * R_st    * R_st  ) &
    + ((dF_dpsi_dZ2 * psi_s) + (dF_dR_dZ2 * R_s) + (dF_dZ3 * Z_s))             * (Z_s    * Z_tt    + 2.0 * Z_st    * Z_t   ) &
    +                           dF_dZ2                                         * (Z_ss   * Z_tt    + 2.0 * Z_sst   * Z_t   ) &
    +                           dF_dZ2                                         * (Z_s    * Z_stt   + 2.0 * Z_st    * Z_st  ) &
    + ((dF_dpsi2_dR * psi_s) + (dF_dpsi_dR2 * R_s) + (dF_dpsi_dR_dZ * Z_s)) &
                                                         * (R_s  * psi_tt  + psi_s  * R_tt  + 2.0 * R_st  * psi_t  + 2.0 * R_t  * psi_st ) &
    +                           dF_dpsi_dR               * (R_ss * psi_tt  + psi_ss * R_tt  + 2.0 * R_sst * psi_t  + 2.0 * R_st * psi_st ) &
    +                           dF_dpsi_dR               * (R_s  * psi_stt + psi_s  * R_stt + 2.0 * R_st  * psi_st + 2.0 * R_t  * psi_sst) &
    + ((dF_dpsi2_dZ * psi_s) + (dF_dpsi_dR_dZ * R_s) + (dF_dpsi_dZ2 * Z_s)) &
                                                         * (Z_s  * psi_tt  + psi_s  * Z_tt  + 2.0 * Z_st  * psi_t  + 2.0 * Z_t  * psi_st ) &
    +                           dF_dpsi_dZ               * (Z_ss * psi_tt  + psi_ss * Z_tt  + 2.0 * Z_sst * psi_t  + 2.0 * Z_st * psi_st ) &
    +                           dF_dpsi_dZ               * (Z_s  * psi_stt + psi_s  * Z_stt + 2.0 * Z_st  * psi_st + 2.0 * Z_t  * psi_sst) &
    + ((dF_dpsi_dR_dZ * psi_s) + (dF_dR2_dZ * R_s) + (dF_dR_dZ2 * Z_s)) &
                                                         * (Z_s  * R_tt    + R_s    * Z_tt  + 2.0 * Z_st  * R_t    + 2.0 * Z_t  * R_st   ) &
    +                           dF_dR_dZ                 * (Z_ss * R_tt    + R_ss   * Z_tt  + 2.0 * Z_sst * R_t    + 2.0 * Z_st * R_st   ) &
    +                           dF_dR_dZ                 * (Z_s  * R_stt   + R_s    * Z_stt + 2.0 * Z_st  * R_st   + 2.0 * Z_t  * R_sst  ) &
    + ((dF_dpsi2 * psi_s) + (dF_dpsi_dR * R_s) + (dF_dpsi_dZ * Z_s))       * psi_stt  &
    +                           dF_dpsi                                    * psi_sstt &
    + ((dF_dpsi_dR * psi_s) + (dF_dR2 * R_s) + (dF_dR_dZ * Z_s))           * R_stt    &
    +                           dF_dR                                      * R_sstt   &
    + ((dF_dpsi_dZ * psi_s) + (dF_dR_dZ * R_s) + (dF_dZ2 * Z_s))           * Z_stt    &
    +                           dF_dZ                                      * Z_sstt   &
    + var_out(5) &
    + var_out(6)


endif

#ifdef fullmhd
if (i_var .eq. 710) then
node_list%node(i_node)%Fprof_eq(1) = var_out(1)
node_list%node(i_node)%Fprof_eq(2) = var_out(2)
node_list%node(i_node)%Fprof_eq(3) = var_out(3)
node_list%node(i_node)%Fprof_eq(4) = var_out(4)
node_list%node(i_node)%Fprof_eq(5) = var_out(5)
node_list%node(i_node)%Fprof_eq(6) = var_out(6)
node_list%node(i_node)%Fprof_eq(7) = var_out(7)
node_list%node(i_node)%Fprof_eq(8) = var_out(8)
node_list%node(i_node)%Fprof_eq(9) = var_out(9)
else
#endif
node_list%node(i_node)%values(1,1,i_var) = var_out(1)
node_list%node(i_node)%values(1,2,i_var) = var_out(2)
node_list%node(i_node)%values(1,3,i_var) = var_out(3)
node_list%node(i_node)%values(1,4,i_var) = var_out(4)
node_list%node(i_node)%values(1,5,i_var) = var_out(5)
node_list%node(i_node)%values(1,6,i_var) = var_out(6)
node_list%node(i_node)%values(1,7,i_var) = var_out(7)
node_list%node(i_node)%values(1,8,i_var) = var_out(8)
node_list%node(i_node)%values(1,9,i_var) = var_out(9)
#ifdef fullmhd
endif
#endif




return
end subroutine project_var_on_node

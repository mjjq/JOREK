subroutine element_matrix_fft

* variable list
* ! --- Initialise
* ! --- Implicit scheme \
$\theta = time\_evol\_theta$ \
$\zeta = time\_evol\_zeta*2*tstep/(tstep + tstep\_prev)$\
* ! --- Determine whether using FFT: use_fft = .True. or .False. \
    ```fortran
    if (use_fft) then
        ! --- Not needed in case of FFT
    else
        ! --- (1st and 2nd derivative of) basis functions in toroidal direction will be used: HZ, HZ_p, HZ_pp
    endif
    ```
* ! --- initialize variables using basis functions, like H, H_s, H_t, H_ss, Hst, Htt, HZ, HZ_p, HZ_pp
    ```fortran
    do i=1,n_vertex_max
    do j=1,n_degrees
        do ms=1, n_gauss
        do mt=1, n_gauss
            x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
            ......
            ! --- Equilibrium psi (n=0 only) for sources
            rho_axisym(ms,mt) = rho_axisym(ms,mt) + nodes(i)%values(1,j,var_rho) * element%size(i,j) * H(i,j,ms,mt)  * HZ(1,1)
            ......
        enddo
        enddo
        do ms=1, n_gauss
        do mt=1, n_gauss
            do k=1,n_var
            do in=1,n_tor
            do mp=1,n_plane
                ! --- store variables
                eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                ......
            enddo
            enddo
            enddo
        enddo
        enddo
    enddo
    enddo
    ```
* ! --- approximate estimate of the element length h_e, needed for Shock capturing stabilzation
* ! --- Sources: current, particle, heat_i/e, V, dV_dpsi, dV_dz
* ! --- Main loops \
    <!-- * $\_s = \partial_s$ 
    * $\_t = \partial_t$ 
    * $\_p = \partial_\varphi$ 
    * $\_R = \partial_R$ 
    * $\_Z = \partial_Z$ -->
    ```fortran
    do i=1,n_vertex_max !< maximum number of corners of an element
    do j=1,n_degrees !< degrees of freedom per variable per node in 2D
        do ms=1, n_gauss
        do mt=1, n_gauss
    
            ! xjac -> Jacobian, i.e., det(J) = R_s * Z_t - R_t * Z_s = x_s * y_t - x_t * y_s
            do mp = 1, n_plane
                ! Calculate AR, AZ, A3 (= \psi = R A_\varphi), 
                    ! UR, UZ, Up (U_\varphi), 
                    ! \rho, Ti, Te of equlibrium, 
                    ! and their derivatives
                ! Calculate pi0, pe0, p0, 
                    ! psi_norm, D_prof, ZKi/e_prof, 
                    ! eta_T, deta_dT, d2eta_dT2, eta_R/Z/p, 
                    ! visco_T, dvisco_dT, 
                    ! ZKi/e_par_T, dZKi/e_par_dT
                ! --- Ion-electron energy transfer, not implemented due to numerical reasons
    ```
    #### $\vec{B}_0\cdot\nabla T_{i/e, 0}$
    * $BgradTi=\vec{B}_0\cdot \nabla T_{i0}$
    * $BgradTe=\vec{B}_0\cdot \nabla T_{e0}$
    #### $\vec{V}_0\cdot\nabla T_{i/e, 0}$
    * $UgradTi=\vec{U}_0\cdot \nabla T_{i0}$
    * $UgradTe=\vec{U}_0\cdot \nabla T_{e0}$
    #### $\vec{B}_0\cdot\nabla \rho_0$
    * $BgradRho=\vec{B}_0\cdot \nabla \rho_0$
    #### $\vec{V}_0\cdot\nabla \rho_0$
    * $UgradRho=\vec{U}_0\cdot \nabla \rho_0$
    #### $\vec{B}_0\cdot\nabla P_{e0}$
    * $BgradPe=\vec{B}_0\cdot \nabla P_{e0}$ 
    #### $\vec{V}_0\cdot\nabla \vec{V}_0$, incomplete, the missed terms of ... will be added when using it.
    * $UgradUR=U_{R0}\frac{\partial U_{R0}}{\partial R} + U_{Z0}\frac{\partial U_{R0}}{\partial Z} + \frac{U_{\varphi 0}}{R}\frac{\partial U_{R0}}{\partial \varphi}+...$
    * $UgradUZ=U_{R0}\frac{\partial U_{Z0}}{\partial R} + U_{Z0}\frac{\partial U_{Z0}}{\partial Z} + \frac{U_{\varphi 0}}{R}\frac{\partial U_{Z0}}{\partial \varphi}$
    * $UgradUp=U_{R0}\frac{\partial U_{\varphi 0}}{\partial R} + U_{Z0}\frac{\partial U_{\varphi 0}}{\partial Z} + \frac{U_{\varphi 0}}{R}\frac{\partial U_{\varphi 0}}{\partial \varphi}+...$
    #### $\nabla\cdot\vec{V}_0$
    * $divU = \frac{\partial U_{R0}}{\partial R} + \frac{U_{R0}}{R} + \frac{\partial U_{Z0}}{\partial Z} + \frac{1}{R}\frac{\partial U_{\varphi 0}}{\partial \varphi}$
    #### $\nabla\cdot(\rho_0\vec{V}_0)$
    * $divRhoU  = \rho_0\nabla\cdot\vec{U}_0 + \vec{U}_0\cdot\nabla\rho$
    #### $\vec{V}_{dia,0}$
    * $VdiaR0=\tau_{IC}\frac{F_0}{R\rho_0B_0^2}(B_{Z0}\frac{\partial p_{i0}}{\partial \varphi} - RB_{\varphi 0}\frac{\partial p_{i0}}{\partial Z})$
    * $VdiaZ0=\tau_{IC}\frac{F_0}{R\rho_0B_0^2}(RB_{\varphi 0}\frac{\partial p_{i0}}{\partial R} - B_{R0}\frac{\partial p_{i0}}{\partial \varphi})$
    * $VdiaP0=\tau_{IC}\frac{F_0}{\rho_0B_0^2}(B_{R0}\frac{\partial p_{i0}}{\partial Z} - B_{Z0}\frac{\partial p_{i0}}{\partial R})$
    #### $\vec{V}_{dia,0}\cdot \nabla \vec{V_0}$, incomplete, the missed terms of ... will be added when using it.
    * $VdiaGradUR  = VdiaR0 * \frac{\partial U_{R0}}{\partial R} + VdiaZ0 * \frac{\partial U_{R0}}{\partial Z} + VdiaP0 * \frac{\partial U_{R0}}{\partial \varphi} / R+...$
    * $VdiaGradUZ  = VdiaR0 * \frac{\partial U_{Z0}}{\partial R} + VdiaZ0 * \frac{\partial U_{Z0}}{\partial Z} + VdiaP0 * \frac{\partial U_{Z0}}{\partial \varphi} / R$
    * $VdiaGradUp  = VdiaR0 * \frac{\partial U_{\varphi 0}}{\partial R} + VdiaZ0 * \frac{\partial U_{\varphi 0}}{\partial Z} + VdiaP0 * \frac{\partial U_{\varphi 0}}{\partial \varphi} / R+...$
    #### Toroidal velocity
    * $Vt0 = V_{source,\varphi 0}$
    * $Vt0\_R = \partial V_{source,\varphi 0}/\partial R$
    * $Vt0\_Z = \{\partial_z V\}_{source,0} + \partial V_{\varphi, 0}/\partial Z$
    ```fortran
                if (NEO) then ! not used in AUG flux pumping cases
                    PneoR = ...
                    PneoZ = ...
                else
                    PneoR = 0
                    PneoZ = 0
                endif
    ```
    #### Calculate $\partial/\partial R(Z,\varphi)$ for several equlibrium variables
    * $psieq\_R=\partial \psi_0/\partial R$
    * $psieq\_Z=\partial \psi_0/\partial Z$ \
......

    ```fortran
                do im=n_tor_start, n_tor_end
                    ! --- test functions (V*)
                    ......
                    ! --- Terms for Integration-by-parts
    ```
    #### Terms for Integration-by-parts: $\vec{X}\cdot \nabla v$
    * $XgradVstar\_\_p = X_{R0} \cdot \frac{\partial v}{\partial R} + X_{Z0} \cdot \frac{\partial v}{\partial Z}$ # poloidal component
    * $XgradVstar\_\_k = X_{\varphi 0} \cdot \frac{1}{R}\frac{\partial v}{\partial \varphi}$ # toroidal component
        * $X$ stands for $\vec{B}_0, \vec{V}_0,\nabla\rho_0,\nabla T_{i0},\nabla T_{e0}$
        * $gradVstar$ stands for $\nabla v$
    <!-- * $BgradVstar\_\_p = B_{R0} \cdot \frac{\partial v}{\partial R} + B_{Z0} \cdot \frac{\partial v}{\partial Z}$
    * $BgradVstar\_\_k = B_{\varphi 0} \cdot \frac{1}{R}\frac{\partial v}{\partial \varphi}$ -->
    <!-- * $UgradVstar\_\_p = V_{R0} \cdot \frac{\partial v}{\partial R} + V_{Z0} \cdot \frac{\partial v}{\partial Z}$
    * $UgradVstar\_\_k = V_{\varphi 0} \cdot \frac{1}{R}\frac{\partial v}{\partial \varphi}$ -->

    #### $-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0$
    * $Qconv\_UR=\{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_R$
    * $Qconv\_UZ=\{...\}_Z$
    * $Qconv\_Up=\{...\}_\varphi$

    #### $\vec{J}_0\times\vec{B}_0$ in the 3 components of momentum equation?
    * $JxB\_UR\_\_p=-B_{R0}\cdot \{\vec{B}_0\cdot \nabla v\}_{poloidal} - B_{\varphi 0}^2v/R+(\partial_R v+ v/R)\cdot B_0^2/2$
    * $JxB\_UR\_\_k=-B_{R0}\cdot \{\vec{B}_0\cdot \nabla v\}_{toroidal}$
    * $JxB\_UZ\_\_p=-B_{Z0}\cdot \{\vec{B}_0\cdot \nabla v\}_{poloidal}+\partial_Z v\cdot B_0^2/2$
    * $JxB\_UZ\_\_k=-B_{Z0}\cdot \{\vec{B}_0\cdot \nabla v\}_{toroidal}$
    * $JxB\_Up\_\_p=-B_{\varphi 0}\cdot \{\vec{B}_0\cdot \nabla v\}_{poloidal}+B_{\varphi 0}B_{R0} v/R$
    * $JxB\_Up\_\_k=-B_{\varphi 0}\cdot \{\vec{B}_0\cdot \nabla v\}_{toroidal}+\partial_\varphi v\cdot B_0^2/2R$

    #### Diamagnetic terms: $\vec{V}_{dia,0}\cdot\nabla v$
    * $VdiaGradVstar\_\_p=\vec{V}_{dia,R0}\partial_R v+\vec{V}_{dia,Z0}\partial_Z v$
    * $VdiaGradVstar\_\_k=\vec{V}_{dia,\varphi 0}\partial_\varphi v/R$

    #### Viscous terms: 
    * $Qvisc\_UR\_\_p = - (\partial_Rv + v/R) \nabla\cdot \vec{V}_0 + \partial_Zv (\partial_R V_{Z0} - \partial_Z V_{R0})$
    * $Qvisc\_UR\_\_k = - \partial_\varphi v/R^2 (\partial_\varphi V_{R0} - V_{\varphi 0} - R\partial_R V_{\varphi 0})$

    * $Qvisc\_UZ\_\_p = - \partial_Zv \nabla\cdot \vec{V}_0 + \partial_Rv (\partial_R V_{Z0} - \partial_Z V_{R0})$
    * $Qvisc\_UZ\_\_k = + \partial_\varphi v/R^2 (R\partial_Z V_{\varphi 0} - \partial_\varphi V_{Z 0})$

    * $Qvisc\_Up\_\_p = - \partial_Z v/R(R\partial_ZV_{\varphi 0} - \partial_\varphi V_{Z0}) + (\partial_Rv + v/R)/R(\partial_\varphi V_{R0} - V_{\varphi 0} - R\partial_RV_{\varphi0})$
    * $Qvisc\_Up\_\_k = - \partial_\varphi v/R \nabla\cdot \vec{V}_0$

    #### The toroidal momentum source:
    * $Qvisc\_UR\_\_k+=\partial_\varphi v/R^2 (- V_{source,\varphi 0} - R\partial_RV_{source,\varphi 0})$
    * $Qvisc\_UZ\_\_k-= \partial_\varphi v/R^2(R\partial_ZV_{source,\varphi 0})$
    * $\color{red}Qvisc\_Up\_\_p+=\partial_Zv/R (R\partial_ZV_{source,\varphi 0}) - (\partial_Rv + v/R)/R (- V_{source,\varphi 0} - R\partial_RV_{source,\varphi 0})$

    #### fourth order diffusion terms: $lap\_X=\Delta^\star X=\partial_RX / R + \partial_R^2X + \partial_Z^2X$, the Laplacian operator $\Delta^\star$ doesn't include the toroidal component
    * where $X$ stands for $v,\ A_{R0},\ A_{Z0},\ A_{30},\ U_{R0},\ U_{Z0},\ U_{\varphi 0},\ \rho_0,\ T_{i0},\ T_{e0}$
    * e.g., $lap\_Vstar = \partial_Rv / R + \partial_R^2v + \partial_Z^2v$

    #### calculate rhs_p_ij, rhs_k_ij, Pvec_prev, Qvec_p, Qvec_k
    * rhs_p_ij  = 0
    * rhs_k_ij  = 0
    * Pvec_prev = 0 ! The time derivative part
    * Qvec_p    = 0 ! The rest of the RHS (poloidal part)
    * Qvec_k    = 0 ! The rest of the RHS (toroidal part that has phi-derivatives of the test-function)
    * delta_g: Change of variable values and derivatives in last timestep, $\color{red}\partial_t$ ?
    * equation 1 (R component induction equation):
        * $Pvec\_prev(A_R)=v\partial_tA_R$
        * $Qvec\_p(A_R)=+v(V_{Z0}B_{\varphi 0}-V_{\varphi 0}B_{Z0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\tau_{IC,A_{RZ}}v\tau_{IC}\frac{F_0}{\rho_0B_0^2}B_{R0}\vec{B}_0\cdot\nabla P_{e0}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\eta_{A_{RZ}}v(\partial_Z\eta B_{\varphi 0}-\partial_\varphi \eta B_{Z0}/R)$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{A_{RZ}}\eta_T(-\partial_ZvB_{\varphi 0}) $\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ -\eta_{A_{RZ,constant}}(-\partial_ZvB_{\varphi 0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\eta_{A_{RZ}}\eta_TvJ_{source,R}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\eta_{A_{RZ}}\eta_{num}\Delta^\star v\Delta^\star A_{R0}$
        * $Qvec\_k(A_R)=-\eta_{A_{RZ}}\eta_T(\partial_\varphi vB_{Z0}/R)$

    * equation 2 (Z component induction equation):
        * $Pvec\_prev(A_Z) =   v \partial_t A_Z $
        * $Qvec\_p(A_Z) = + v(V_{\varphi 0}B_{R 0}-V_{R 0}B_{\varphi 0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\tau_{IC,A_{RZ}}v\tau_{IC}\frac{F_0}{\rho_0B_0^2}B_{Z0}\vec{B}_0\cdot\nabla P_{e0}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\eta_{A_{RZ}}v(\partial_\varphi\eta B_{R0}/R-\partial_R \eta B_{\varphi 0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +-\eta_{A_{RZ}}\eta_T(\partial_RvB_{\varphi 0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ -\eta_{A_{RZ,constant}}(\partial_RvB_{\varphi 0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\eta_{A_{RZ}}\eta_TvJ_{source,Z}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\eta_{A_{RZ}}\eta_{num}\Delta^\star v\Delta^\star A_{Z0}$
        * $Qvec\_k(A_Z)=-\eta_{A_{RZ}}\eta_T(-\partial_\varphi vB_{R0}/R)$

    * equation 3 (PHI component induction equation):
        * $Pvec\_prev(A_3) =   v \partial_t A_3 $
        * $Qvec\_p(A_3) = + Rv(V_{R 0}B_{Z 0}-V_{Z 0}B_{R 0})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\tau_{IC,A_{\varphi}}v\tau_{IC}\frac{F_0}{\rho_0B_0^2}RB_{\varphi 0}\vec{B}_0\cdot\nabla P_{e0}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ -\eta_T\partial_ZvRB_{R0}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \eta_T(2v+R\partial_Rv)B_{Z0}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + Rv(\partial_R\eta B_{Z0}-\partial_Z\eta B_{R0})$ \
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \eta_T v J_{source,\varphi}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ -\eta_{num}\Delta^\star v\Delta^\star A_{30}$

    * equation 4 (R component momentum equation):
        * $Pvec\_prev(V_R) =   v\rho_0 \partial_tV_R$
        * $Qvec\_p(V_R)=+\{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_R$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +(\partial_R v+v/R)P_0$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \underline{JxB\_UR\_\_p}$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \mu_T\underline{Qvisc\_UR\_\_p}$ \
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - v\sout{PneoR}$ \
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - v S_{particle}V_{R0}$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - \mu_{num}\Delta^\star v\Delta^\star V_{R0}$
        * $Qvec\_k(V_R)=+\underline{JxB\_UR\_\_k}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \mu_T \underline{Qvisc\_UR\_\_k}$

    * equation 5 (Z component momentum equation):
        * $Pvec\_prev(V_Z) =   v\rho_0 \partial_tV_Z$
        * $Qvec\_p(V_Z)=+\{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_Z$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +\partial_Z vP_0$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \underline{JxB\_UZ\_\_p}$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \mu_T\underline{Qvisc\_UZ\_\_p}$ \
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - v\sout{PneoZ}$ \
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - v S_{particle}V_{Z0}$\
         $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - \mu_{num}\Delta^\star v\Delta^\star V_{Z0}$
        * $Qvec\_k(V_Z)=+\underline{JxB\_UZ\_\_k}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + \mu_T \underline{Qvisc\_UZ\_\_k}$

    * equation 6 (Phi component momentum equation):
        * $Pvec\_prev(V_\varphi)=+v\rho_0B_{R0}\partial_tV_R$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +v\rho_0B_{Z0}\partial_tV_Z$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +v\rho_0B_{\varphi 0}\partial_tV_\varphi$
        * $Qvec\_p(V_\varphi)=+\vec{B}_0\cdot\{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + B_{R0}\mu_T\underline{Qvisc\_UR\_\_p}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + B_{Z0}\mu_T\underline{Qvisc\_UZ\_\_p}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + B_{\varphi 0}\mu_T\underline{Qvisc\_Up\_\_p}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ + P0\cdot\{\vec{B}_0\cdot\nabla v\}_{poloidal} $\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ - v(\sout{B_{R0}*PneoR + B_{Z0}*PneoZ})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ -vS_{particle}\vec{B}_0\cdot\vec{V}_0$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ -\mu_{num}\Delta^\star v\vec{B}_0\cdot\Delta^\star\vec{V}_0$
        * $Qvec\_k(V_\varphi)=+B_{R0}\mu_T\underline{Qvisc\_UR\_k}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +B_{Z0}\mu_T\underline{Qvisc\_UZ\_k}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +B_{\varphi 0}\mu_T\underline{Qvisc\_Up\_k}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ +P_{0}\{\vec{B}_0\cdot\nabla v\}_{toroidal}$

    * equation 7 (Density equation):
        * $Pvec\_prev(\rho)=v\partial_t\rho$
        * $Qvec\_p(\rho)=-v(\rho_0\nabla\cdot \vec{V}_0+\vec{V}\cdot\nabla\rho_0)$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ +\rho_0\{\vec{V}_{dia,0}\cdot\nabla v\}_{poloidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -D\{\nabla\rho_0\cdot\nabla v\}_{poloidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -(D_\parallel-D)\{\vec{B}_0\cdot\nabla v\}_{poloidal}\vec{B}_0\cdot\nabla\rho_0/B_0^2$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ +vS_{particle}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -D_{\perp, num}\Delta^\star v\Delta^\star\rho_0$
        * $Qvec\_k(\rho)=+\rho_0\{\vec{V}_{dia,0}\cdot\nabla v\}_{toroidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -D\{\nabla\rho_0\cdot\nabla v\}_{toroidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -(D_\parallel-D)\{\vec{B}_0\cdot\nabla v\}_{toroidal}\vec{B}_0\cdot\nabla \rho_0/B_0^2$
    * equation 8 (Ion Pressure equation):
        * $Pvec\_prev(T_i)=\ \ v\rho_0\partial_tT_i$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +vT_{i0}\partial_t\rho$
        * $Qvec\_p(T_i) = +v(-\rho_0\vec{V}_0\cdot\nabla T_{i0}-T_{i0}\vec{V}_0\cdot\nabla\rho_0-\gamma P_{i0}\nabla\cdot\vec{V}_0)$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ +vS_{heat,i}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ +v(\gamma-1)\sout{Qvisc\_T} $\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -\kappa_i\{\nabla T_{i0}\cdot\nabla v\}_{poloidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -(\kappa_{i,T,\parallel}-\kappa_i)\{\vec{B}_0\cdot\nabla v\}_{poloidal}\vec{B}_0\cdot\nabla T_{i0}/B_0^2 $\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ + v\delta T_{i,e}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ - \kappa_{i,num,\perp} \Delta^\star v\Delta^\star T_{i0}$
        * $Qvec\_k(T_i)=-\kappa_i\{\nabla T_{i0}\cdot\nabla v\}_{toroidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -(\kappa_{i,T,\parallel}-\kappa_i)\{\vec{B}_0\cdot\nabla v\}_{toroidal} \vec{B}_0\cdot\nabla T_{i0}/B_0^2$
    * equation 9 (Electron Pressure equation):
        * $Pvec\_prev(T_e)=\ \ v\rho_0\partial_tT_e$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +vT_{e0}\partial_t\rho$
        * $Qvec\_p(T_e) = +v(-\rho_0\vec{V}_0\cdot\nabla T_{e0}-T_{e0}\vec{V}_0\cdot\nabla\rho_0-\gamma P_{e0}\nabla\cdot\vec{V}_0)$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ +vS_{heat,e}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ +v(\gamma-1)\sout{Qvisc\_T} $\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -\kappa_e\{\nabla T_{e0}\cdot\nabla v\}_{poloidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -(\kappa_{e,T,\parallel}-\kappa_e)\{\vec{B}_0\cdot\nabla v\}_{poloidal}\vec{B}_0\cdot\nabla T_{e0}/B_0^2 $\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ + v\delta T_{e,i}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ - \kappa_{e,num,\perp} \Delta^\star v\Delta^\star T_{e0}$
        * $Qvec\_k(T_e)=-\kappa_i\{\nabla T_{e0}\cdot\nabla v\}_{toroidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ -(\kappa_{e,T,\parallel}-\kappa_e)\{\vec{B}_0\cdot\nabla v\}_{toroidal} \vec{B}_0\cdot\nabla T_{e0}/B_0^2$
    
    ```fortran
                    ! --- Fill Up the RHS
                    if (use_fft) then
                    index_ij =       n_var*n_degrees*(i-1) +       n_var*(j-1) + 1
                    do ivar= 1,n_var
                        RHS_p_ij(ivar) = tstep * Qvec_p(ivar) + zeta * Pvec_prev(ivar) * tstep / tstep_prev
                        RHS_k_ij(ivar) = tstep * Qvec_k(ivar)

                        ij = index_ij + (ivar-1)
                        RHS_p(mp,ij)   =  RHS_p(mp,ij) + RHS_p_ij(ivar) * wst * R * xjac
                        RHS_k(mp,ij)   =  RHS_k(mp,ij) + RHS_k_ij(ivar) * wst * R * xjac
                    enddo
                    else
                    index_ij = n_tor_local*n_var*n_degrees*(i-1) + n_tor_local*n_var*(j-1) + im - n_tor_start +1
                    do ivar= 1,n_var
                        RHS_p_ij(ivar) = tstep * Qvec_p(ivar) + zeta * Pvec_prev(ivar) * tstep / tstep_prev
                        RHS_k_ij(ivar) = tstep * Qvec_k(ivar)

                        ij = index_ij + (ivar-1)*n_tor_local
                        RHS(ij)        =  RHS(ij) + (RHS_p_ij(ivar) + RHS_k_ij(ivar)) * wst * R * xjac
                    enddo
                    endif
    ```

    ```fortran
                    do k=1,n_vertex_max
                    do l=1,n_degrees
                    do in =  n_tor_start, n_tor_end
                    ! --- Basis functions
                    ......
                    ! --- Linearised quantities
                    ......
                    amat      = 0
                    Pjac      = 0 ! time derivative part
                    Qjac_p    = 0 ! rest of the LHS (poloidal part)
                    Qjac_k    = 0 ! rest of the LHS (toroidal part with phi-derivatives of the test-function)
                    Qjac_n    = 0 ! rest of the LHS (toroidal part with phi-derivatives of the basis-functions)
                    Qjac_kn   = 0 ! rest of the LHS (toroidal part with phi-derivatives of the basis-functions and the test-functions)
                    Qjac_pnn  = 0
    ```
    * Hint: why $\partial B_\varphi/\partial A_R = -\partial_Z v$, i.e., in the code: Bp0_AR = -AR_Z = -bf_Z?
        * $B_\varphi=-\partial_ZA_R+\partial_RA_Z$
        * $\dfrac{\partial(vB_\varphi)}{\partial A_R}=\dfrac{\partial(-v\partial_Z A_R)}{\partial A_R}=-\dfrac{\partial v}{\partial A_R}\dfrac{\partial A_R}{\partial Z}-v\dfrac{\partial^2A_R}{\partial Z\partial A_R}=-\dfrac{\partial v}{\partial Z}$
        * $\frac{\partial[\Delta^\star (v A_{R0})]}{\partial A_R}=\Delta^*v=\Delta^\star b_f=lap\_bf$
        * Similiar for other quantities.
        * The equations below are nonlinear, e.g., $v(V_{Z0} * \frac{\partial B_{\varphi 0}}{\partial A_R} )$ contains the product of $v$ and $\partial_Zv$
    * equation 1 (R component induction equation):
        * $Pjac   (A_R, A_R) = v * A_R$
        * $Qjac_p (A_R,A_R) = + v(V_{Z0} \frac{\partial B_{\varphi 0}}{\partial A_R} )$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +v\tau_{IC}\frac{F_0}{\rho_0B_0^2}B_{R0}\{\frac{\partial (\vec{B}_0\cdot \nabla P_e)}{\partial A_R}\}_{poloidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -v\tau_{IC}\frac{F_0}{\rho_0B_0^4}B_{R0}\vec{B}_0\cdot\nabla P_e\{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +v(\partial_Z\eta\frac{\partial B_{\varphi 0}}{\partial A_R})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -\eta_T(-\partial_Zv\frac{\partial B_{\varphi 0}}{\partial A_R})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -\eta_{constant}(-\partial_Zv\frac{\partial B_{\varphi 0}}{\partial A_R})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -\eta_{num}\Delta^\star v\Delta^\star b_f$
        * $Qjac\_n (A_R, A_R)=+v(-V_{\varphi 0}\frac{\partial B_{Z0}}{\partial A_R})$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +v\tau_{IC}\frac{F_0}{\rho_0B_0^2}B_{R0}\{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_R}\}_{toroidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -v\tau_{IC}\frac{F_0}{\rho_0B_0^4}B_{R0}\vec{B}_0\cdot\nabla P_e\{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal}$\
        $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +v(-\partial_\varphi\eta\frac{\partial B_{Z0}}{\partial A_R}\frac{1}{R})$\
        * $Qjac\_kn(A_R,A_R)=-  \eta_T * ( + {\partial_\varphi v} \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} / R)$
        * $  Qjac\_p (A_R,A_Z) = + v * (V_{Z0} * \frac{\partial B_{\varphi 0}}{\partial A_Z} - V_{\varphi 0} * \frac{\partial B_{Z0}}{\partial A_Z})          $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{R0} * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{R0} * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\partial_Z\eta * \frac{\partial B_{\varphi 0}}{\partial A_Z} - \partial_\varphi\eta * \frac{\partial B_{Z0}}{\partial A_Z} / R ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \eta_T * ( - \partial_Z v * \frac{\partial B_{\varphi 0}}{\partial A_Z} ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{constant} ( - \partial_Z v * \frac{\partial B_{\varphi 0}}{\partial A_Z} ) $
        * $  Qjac\_n (A_R,A_Z) = +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} * (\vec{B}_0\cdot\nabla P_e) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{R0}       * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{R0}       * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $
        * $  Qjac\_k (A_R,A_Z) = -  \eta_T * ( + {\partial_\varphi v} * \frac{\partial B_{Z0}}{\partial A_Z} / R) $
        * $  Qjac\_p (A_R,A_3) = + v * (V_{Z0} * \frac{\partial B_{\varphi 0}}{\partial A_3} - V_{\varphi 0} * \frac{\partial B_{Z 0}}{\partial A_3})          $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * \frac{\partial B_{R 0}}{\partial A_3} * (\vec{B}_0\cdot\nabla P_e) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{R0}    * \frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{R0}    * (\vec{B}_0\cdot\nabla P_e) * \frac{\partial B_0^2}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\partial_Z\eta * \frac{\partial B_{\varphi 0}}{\partial A_3} - \partial_\varphi\eta * \frac{\partial B_{Z 0}}{\partial A_3} / R ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \eta_T * ( - \partial_Z v * \frac{\partial B_{\varphi 0}}{\partial A_3} ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{constant} ( - \partial_Z v * \frac{\partial B_{\varphi 0}}{\partial A_3} ) $
        * $  Qjac\_k (A_R,A_3) = -  \eta_T * ( + {\partial_\varphi v} * \frac{\partial B_{Z 0}}{\partial A_3} / R)$
        * $  Qjac\_p (A_R,V_Z) = + v * (  V_Z * B_{\varphi 0}) $
        * $  Qjac\_p (A_R,V_\varphi) = + v * (- V_\varphi * B_{Z0})$
        * $  Qjac\_p (A_R,\rho)= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{R0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial\rho}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0^2B_0^2} * B_{R0} * (\vec{B}_0\cdot\nabla P_e) * \rho$
        * $  Qjac\_n (A_R,\rho)= -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{R0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial \rho}\}_{toroidal} $
        * $  Qjac\_p (A_R,T_e )= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{R0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial T_e}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\frac{\partial}{\partial Z}(\frac{\partial \eta}{\partial T} T_e) * B_{\varphi 00} - \{\frac{\partial}{\partial \varphi}(\frac{\partial \eta}{\partial T} T_e)\}_{poloidal} * B_{Z0} / R ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \frac{\partial \eta}{\partial T} T_e * ( - \partial_Z v * B_{\varphi 00} )                   $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  \frac{\partial \eta}{\partial T} T_e * v * S_{JR} $
        * $  Qjac\_n (A_R,T_e )= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{R0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial T_e}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (              - \{\frac{\partial}{\partial \varphi}(\frac{\partial \eta}{\partial T} T_e)\}_{toroidal} * B_{Z0} / R ) $
        * $  Qjac\_k (A_R,T_e )= -  \frac{\partial \eta}{\partial T} T_e * ( + {\partial_\varphi v} * B_{Z0} / R) $
    * equation 2   (Z component induction equation)
        * $  Pjac   (A_Z,A_Z) =   v * A_Z$
        * $  Qjac\_p (A_Z,A_R) = + v * (V_{\varphi 0} * \frac{\partial B_{R0}}{\partial A_R} - V_{R0} * \frac{\partial B_{\varphi 0}}{\partial A_R})         $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{Z0} * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_R}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{Z0} * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\partial_\varphi\eta / R * \frac{\partial B_{R0}}{\partial A_R} - \partial_R \eta * \frac{\partial B_{\varphi 0}}{\partial A_R}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \eta_T * ( + v_R * \frac{\partial B_{\varphi 0}}{\partial A_R} ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{constant} ( + v_R * \frac{\partial B_{\varphi 0}}{\partial A_R} )$
        * $  Qjac\_n (A_Z,A_R) = +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} * (\vec{B}_0\cdot\nabla P_e) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{Z0}       * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_R}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{Z0}       * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $
        * $  Qjac\_k (A_Z,A_R) = -  \eta_T * ( - {\partial_\varphi v} * \frac{\partial B_{R0}}{\partial A_R} / R) $
        * $  Qjac\_p (A_Z,A_Z) = + v * (- V_{R0}   * \frac{\partial B_{\varphi 0}}{\partial A_Z}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{Z0} * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{Z0} * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (- \partial_R \eta * \frac{\partial B_{\varphi 0}}{\partial A_Z}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \eta_T * ( + v_R * \frac{\partial B_{\varphi 0}}{\partial A_Z} ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{constant} ( + v_R * \frac{\partial B_{\varphi 0}}{\partial A_Z} ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \eta_{num} * \Delta^\star v * \Delta^\star b_f $
        * $  Qjac\_n (A_Z,A_Z) = + v * (V_{\varphi 0} * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal})       $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{Z0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{Z0} * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\partial_\varphi\eta / R * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal}) $
        * $  Qjac\_kn(A_Z,A_Z) = -  \eta_T * ( - {\partial_\varphi v} * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} / R) $
        * $  Qjac\_p (A_Z,A_3) = + v * (V_{\varphi 0} * \frac{\partial B_{R 0}}{\partial A_3} - V_{R0} * \frac{\partial B_{\varphi 0}}{\partial A_3})         $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * \frac{\partial B_{Z 0}}{\partial A_3} * (\vec{B}_0\cdot\nabla P_e) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * B_{Z0}    * \frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * B_{Z0}    * (\vec{B}_0\cdot\nabla P_e) * \frac{\partial B_0^2}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\partial_\varphi\eta / R * \frac{\partial B_{R 0}}{\partial A_3} - \partial_R \eta * \frac{\partial B_{\varphi 0}}{\partial A_3}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \eta_T * ( + v_R * \frac{\partial B_{\varphi 0}}{\partial A_3} ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{constant} ( + v_R * \frac{\partial B_{\varphi 0}}{\partial A_3} ) $
        * $  Qjac\_k (A_Z,A_3) = -  \eta_T * ( - {\partial_\varphi v} * \frac{\partial B_{R 0}}{\partial A_3} / R) $
        * $  Qjac\_p (A_Z,V_R) = + v * (- V_R * B_{\varphi 0}) $
        * $  Qjac\_p (A_Z,V_\varphi) = + v * (  V_\varphi * B_{R0}) $
        * $  Qjac\_p (A_Z,\rho)= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{Z0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial\rho}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0^2B_0^2} * B_{Z0} * (\vec{B}_0\cdot\nabla P_e) * \rho $
        * $  Qjac\_n (A_Z,\rho)= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{Z0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial \rho}\}_{toroidal} $
        * $  Qjac\_p (A_Z,T_e )= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{Z0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial T_e}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\{\frac{\partial}{\partial \varphi}(\frac{\partial \eta}{\partial T} T_e)\}_{poloidal} / R * B_{R0} - \frac{\partial}{\partial R}(\frac{\partial \eta}{\partial T} T_e) * B_{\varphi 00}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  \frac{\partial \eta}{\partial T} T_e * ( + v_R * B_{\varphi 00} )                  $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  \frac{\partial \eta}{\partial T} T_e * v * S_{JZ} $
        * $  Qjac\_n (A_Z,T_e )= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * B_{Z0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial T_e}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * (\{\frac{\partial}{\partial \varphi}(\frac{\partial \eta}{\partial T} T_e)\}_{toroidal} / R * B_{R0} ) $
        * $  Qjac\_k (A_Z,T_e )= -  \frac{\partial \eta}{\partial T} T_e * ( - {\partial_\varphi v} * B_{R0} / R) $
    * equation 3   (Phi component induction equation)
        * $  Pjac   (A_3,A_3) =   v * A3 $
        * $  Qjac\_p (A_3,A_R) = + R * v * (- V_{Z0} * \frac{\partial B_{R0}}{\partial A_R})                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*\frac{\partial B_{\varphi 0}}{\partial A_R} * (\vec{B}_0\cdot\nabla P_e) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*B_{\varphi 0}    * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_R}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * R*B_{\varphi 0}    * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_T * \partial_Z v * R                * \frac{\partial B_{R0}}{\partial A_R} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + R * v * (- \partial_Z\eta * \frac{\partial B_{R0}}{\partial A_R}) $
        * $  Qjac\_n (A_3,A_R) = + R * v * (V_{R0} * \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal})                  $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*B_{\varphi 0}    * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_R}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * R*B_{\varphi 0}    * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + \eta_T * ( 2 * v + R * v_R ) * \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + R * v * (\partial_R \eta * \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal}) $
        * $  Qjac\_p (A_3,A_Z) = + R * v * (V_{R0} * \frac{\partial B_{Z0}}{\partial A_Z})                  $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*\frac{\partial B_{\varphi 0}}{\partial A_Z} * (\vec{B}_0\cdot\nabla P_e) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*B_{\varphi 0}    * \{\frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * R*B_{\varphi 0}    * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + \eta_T * ( 2 * v + R * v_R ) * \frac{\partial B_{Z0}}{\partial A_Z} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + R * v * (\partial_R \eta * \frac{\partial B_{Z0}}{\partial A_Z}) $
        * $  Qjac\_n (A_3,A_Z) = + R * v * (- V_{Z0} * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal})                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*B_{\varphi 0}    * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * R*B_{\varphi 0}    * (\vec{B}_0\cdot\nabla P_e) * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_T * \partial_Z v * R                * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + R * v * (- \partial_Z\eta * \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal}) $
        * $  Qjac\_p (A_3,A_3) = + R * v * (V_{R0} * \frac{\partial B_{Z 0}}{\partial A_3} - V_{Z0} * \frac{\partial B_{R 0}}{\partial A_3})   $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2}    * R*B_{\varphi 0}    * \frac{\partial (\vec{B}_0\cdot\nabla P_e)}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^4} * R*B_{\varphi 0}    * (\vec{B}_0\cdot\nabla P_e) * \frac{\partial B_0^2}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_T * \partial_Z v * R                * \frac{\partial B_{R 0}}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + \eta_T * ( 2 * v + R * v_R ) * \frac{\partial B_{Z 0}}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + R * v * (\partial_R \eta * \frac{\partial B_{Z 0}}{\partial A_3} - \partial_Z\eta * \frac{\partial B_{R 0}}{\partial A_3}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \eta_{num} * \Delta^\star v * \Delta^\star b_f $
        * $  Qjac\_p (A_3,V_R) = + R * v * (  V_R * B_{Z0}) $
        * $  Qjac\_p (A_3,V_Z) = + R * v * (- V_Z * B_{R0}) $
        * $  Qjac\_p (A_3,\rho)= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * R*B_{\varphi 0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial\rho}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -  v * \tau_{IC}*\frac{F_0}{\rho_0^2B_0^2} * R*B_{\varphi 0} * (\vec{B}_0\cdot\nabla P_e) * \rho $
        * $  Qjac\_n (A_3,\rho)= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * R*B_{\varphi 0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial \rho}\}_{toroidal} $
        * $  Qjac\_p (A_3,T_e )= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * R*B_{\varphi 0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial T_e}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - \frac{\partial \eta}{\partial T} T_e * \partial_Z v * R                * B_{R0} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + \frac{\partial \eta}{\partial T} T_e * ( 2 * v + R * v_R ) * B_{Z0} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + \frac{\partial \eta}{\partial T} T_e * v * S_{J\varphi} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ + R * v * (\frac{\partial}{\partial R}(\frac{\partial \eta}{\partial T} T_e) * B_{Z0} - \frac{\partial}{\partial Z}(\frac{\partial \eta}{\partial T} T_e) * B_{R0}) $
        * $  Qjac\_n (A_3,T_e )= +  v * \tau_{IC}*\frac{F_0}{\rho_0B_0^2} * R*B_{\varphi 0} * \{\frac{\partial(\vec{B}_0\cdot\nabla P_e)}{\partial T_e}\}_{toroidal} $
    * equation 4   (R component momentum equation)

        * $  Pjac(V_R,V_R )    =   v * \rho_0 * V_R $
				
        * $  Qjac\_p (V_R,A_R ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_R}\}_{poloidal} + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_R}\}_{poloidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial A_R}\}_{poloidal} $

        * $  Qjac\_n (V_R,A_R ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_R}\}_{toroidal} + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_R}\}_{toroidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial A_R}\}_{toroidal} $

        * $  Qjac\_k (V_R,A_R ) =                 + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_R}\}_{toroidal, k} $

        * $  Qjac\_kn(V_R,A_R ) =                 + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_R}\}_{toroidal, kn} $
				
        * $  Qjac\_p (V_R,A_Z ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_Z}\}_{poloidal} + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_Z}\}_{poloidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial A_Z}\}_{poloidal} $
                
        * $  Qjac\_n (V_R,A_Z ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_Z}\}_{toroidal} + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_Z}\}_{toroidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial A_Z}\}_{toroidal} $
                
        * $  Qjac\_k (V_R,A_Z ) =                 + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_Z}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_R,A_Z ) =                 + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_Z}\}_{toroidal, kn} $
				
        * $  Qjac\_p (V_R,A_3 ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_3}\}_{poloidal} + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_3}\}_{poloidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial A_3}\}_{poloidal} $
                
        * $  Qjac\_n (V_R,A_3 ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_3}\}_{toroidal} + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_3}\}_{toroidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial A_3}\}_{toroidal} $
                
        * $  Qjac\_k (V_R,A_3 ) =                 + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_3}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_R,A_3 ) =                 + \{\frac{\partial(\vec{J}\times\vec{B})_R}{\partial A_3}\}_{toroidal, kn} $
				
        * $  Qjac\_p (V_R,V_R ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_R}\}_{poloidal} + \mu_T * Qvisc\_UR\_UR_{poloidal}  - v * \{\frac{\partial \sout{PneoR}}{\partial V_R}\} - v * S_p * V_R $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * \Delta^\star b_f $

        * $  Qjac\_n (V_R,V_R ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_R}\}_{toroidal} + \mu_T * Qvisc\_UR\_UR_{toroidal} $
                
        * $  Qjac\_k (V_R,V_R ) =                 + \mu_T * Qvisc\_UR\_UR_{toroidal, k} $
                
        * $  Qjac\_kn(V_R,V_R ) =                 + \mu_T * Qvisc\_UR\_UR_{toroidal, kn} $
				
        * $  Qjac\_p (V_R,V_Z ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_Z}\}_{poloidal} + \mu_T * Qvisc\_UR\_UR_{poloidal}  - v * $\{\frac{\partial \sout{PneoR}}{\partial V_Z}\}$ $
                
        * $  Qjac\_n (V_R,V_Z ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_Z}\}_{toroidal} + \mu_T * Qvisc\_UR\_UR_{toroidal} $
                
        * $  Qjac\_k (V_R,V_Z ) =                 + \mu_T * Qvisc\_UR\_UR_{toroidal, k} $
                
        * $  Qjac\_kn(V_R,V_Z ) =                 + \mu_T * Qvisc\_UR\_UR_{toroidal, kn} $
				
        * $  Qjac\_p (V_R,V_\varphi ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_\varphi}\}_{poloidal} + \mu_T * Qvisc\_UR\_Up_{poloidal} $
                
        * $  Qjac\_n (V_R,V_\varphi ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_\varphi}\}_{toroidal} + \mu_T * Qvisc\_UR\_Up_{toroidal} $
                
        * $  Qjac\_k (V_R,V_\varphi ) =                 + \mu_T * Qvisc\_UR\_Up_{toroidal, k} $
                
        * $  Qjac\_kn(V_R,V_\varphi ) =                 + \mu_T * Qvisc\_UR\_Up_{toroidal, kn} $
				
        * $  Qjac\_p (V_R,\rho) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial \rho}\}_{poloidal} + ( v_R + v / R ) * (\rho*(T_{i0}+T_{e0})) - v * \{\frac{\partial \sout{PneoR}}{\partial \rho}\}_{poloidal} $
                
        * $  Qjac\_n (V_R,\rho) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial \rho}\}_{toroidal}                                     - v * \{\frac{\partial \sout{PneoR}}{\partial \rho}\}_{toroidal} $
				
        * $  Qjac\_p (V_R,T_i ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial T_i}\}_{poloidal}  + ( v_R + v / R ) * (\rho_0*Ti      ) - v * \{\frac{\partial \sout{PneoR}}{\partial T_i}\}_{poloidal} $
                
        * $  Qjac\_n (V_R,T_i ) =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial T_i}\}_{toroidal}                                      - v * \{\frac{\partial \sout{PneoR}}{\partial T_i}\}_{toroidal} $
				
        * $  Qjac\_p (V_R,T_e ) =                  + ( v_R + v / R ) * (\rho_0*T_e      ) + \partial_T \mu * T_e * Qvisc\_UR_{poloidal} $
                
        * $  Qjac\_k (V_R,T_e ) =                                                      + \partial_T \mu * T_e * Qvisc\_UR_{toroidal, k} $
	
    * equation 5   (Z component momentum equation)

        * $  Pjac(V_Z,V_Z)      =   v * \rho_0 * V_Z $
				
        * $  Qjac\_p (V_Z,A_R )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_R}\}_{poloidal} + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_R}\}_{poloidal}  - v * \{\frac{\partial \sout{PneoZ}}{\partial A_R}\}_{poloidal} $
               
        * $  Qjac\_n (V_Z,A_R )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_R}\}_{toroidal} + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_R}\}_{toroidal}  - v * \{\frac{\partial \sout{PneoZ}}{\partial A_R}\}_{toroidal} $
                
        * $  Qjac\_k (V_Z,A_R )  =                 + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_R}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_Z,A_R )  =                 + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_R}\}_{toroidal, kn} $
				
        * $  Qjac\_p (V_Z,A_Z )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_Z}\}_{poloidal} + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_Z}\}_{poloidal}  - v * \{\frac{\partial \sout{PneoZ}}{\partial A_Z}\}_{poloidal} $
                
        * $  Qjac\_n (V_Z,A_Z )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_Z}\}_{toroidal} + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_Z}\}_{toroidal}  - v * \{\frac{\partial \sout{PneoZ}}{\partial A_Z}\}_{toroidal} $
                
        * $  Qjac\_k (V_Z,A_Z )  =                 + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_Z}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_Z,A_Z )  =                 + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_Z}\}_{toroidal, kn} $
				
        * $  Qjac\_p (V_Z,A_3 )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_3}\}_{poloidal} + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_3}\}_{poloidal}  - v * \{\frac{\partial \sout{PneoZ}}{\partial A_3}\}_{poloidal} $

        * $  Qjac\_n (V_Z,A_3 )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_3}\}_{toroidal} + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_3}\}_{toroidal}  - v * \{\frac{\partial \sout{PneoZ}}{\partial A_3}\}_{toroidal} $
                
        * $  Qjac\_k (V_Z,A_3 )  =                 + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_3}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_Z,A_3 )  =                 + \{\frac{\partial(\vec{J}\times\vec{B})_Z}{\partial A_3}\}_{toroidal, kn} $
				
        * $  Qjac\_p (V_Z,V_R )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_R}\}_{poloidal} + \mu_T * Qvisc\_UZ\_UR_{poloidal} - v * \{\frac{\partial \sout{PneoZ}}{\partial V_R}\} $
                
        * $  Qjac\_n (V_Z,V_R )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_R}\}_{toroidal} + \mu_T * Qvisc\_UZ\_UR_{toroidal} $
                
        * $  Qjac\_k (V_Z,V_R )  =                 + \mu_T * Qvisc\_UZ\_UR_{toroidal, k} $
                
        * $  Qjac\_kn(V_Z,V_R )  =                 + \mu_T * Qvisc\_UZ\_UR_{toroidal, kn} $

        * $  Qjac\_p (V_Z,V_Z )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_Z}\}_{poloidal} + \mu_T * Qvisc\_UZ\_UZ_{poloidal} - v * \{\frac{\partial \sout{PneoZ}}{\partial V_Z}\} - v * S_p * V_Z $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * \Delta^\star b_f $
               
        * $  Qjac\_n (V_Z,V_Z )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_Z}\}_{toroidal} + \mu_T * Qvisc\_UZ\_UZ_{toroidal} $
                
        * $  Qjac\_k (V_Z,V_Z )  =                 + \mu_T * Qvisc\_UZ\_UZ_{toroidal, k} $
               
        * $  Qjac\_kn(V_Z,V_Z )  =                 + \mu_T * Qvisc\_UZ\_UZ_{toroidal, kn} $
				
        * $  Qjac\_p (V_Z,V_\varphi )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_\varphi}\}_{poloidal} + \mu_T * Qvisc\_UZ\_Up_{poloidal} $
                
        * $  Qjac\_n (V_Z,V_\varphi )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_\varphi}\}_{toroidal} + \mu_T * Qvisc\_UZ\_Up_{toroidal} $
                
        * $  Qjac\_k (V_Z,V_\varphi )  =                 + \mu_T * Qvisc\_UZ\_Up_{toroidal, k} $
                
        * $  Qjac\_kn(V_Z,V_\varphi )  =                 + \mu_T * Qvisc\_UZ\_Up_{toroidal, kn} $
				
        * $  Qjac\_p (V_Z,\rho)  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial \rho}\}_{poloidal} + \partial_Z v * (\rho*(T_{i0}+T_{e0})) - v * \{\frac{\partial \sout{PneoZ}}{\partial \rho}\}_{poloidal} $

        * $  Qjac\_n (V_Z,\rho)  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial \rho}\}_{toroidal}                         - v * \{\frac{\partial \sout{PneoZ}}{\partial \rho}\}_{toroidal} $
				
        * $  Qjac\_p (V_Z,T_i )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial T_i}\}_{poloidal}  + \partial_Z v * (\rho_0*Ti      ) - v * \{\frac{\partial \sout{PneoZ}}{\partial T_i}\}_{poloidal} $

        * $  Qjac\_n (V_Z,T_i )  =  -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial T_i}\}_{toroidal}                          - v * \{\frac{\partial \sout{PneoZ}}{\partial T_i}\}_{toroidal} $
				
        * $  Qjac\_p (V_Z,T_e )  =                  + \partial_Z v * (\rho_0*T_e      ) + \partial_T \mu * T_e * Qvisc\_UZ_{poloidal} $

        * $  Qjac\_k (V_Z,T_e )  =                                          + \partial_T \mu * T_e * Qvisc\_UZ_{toroidal, k} $

    * equation 6   (Phi component momentum equation)

        * $  Pjac(V_\varphi,V_R)    = v * B_{R0} * \rho_0 * V_R $

        * $  Pjac(V_\varphi,V_Z)    = v * B_{Z0} * \rho_0 * V_Z $

        * $  Pjac(V_\varphi,V_\varphi)    = v * B_{\varphi 0} * \rho_0 * V_\varphi  $
				
        * $  Qjac\_p (V_\varphi,A_R)  = + \frac{\partial B_{R0}}{\partial A_R} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_R              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_R} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_\varphi              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_R}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_R}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial A_R}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{R0}}{\partial A_R} * \mu_T * Qvisc\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_R} * \mu_T * Qvisc\_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{poloidal}          $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial A_R}\}_{poloidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial A_R}\}_{poloidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (\frac{\partial B_{\varphi 0}}{\partial A_R}*V_{\varphi 0}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * (\frac{\partial B_{R0}}{\partial A_R} * \Delta^\star V_{R0}               + \frac{\partial B_{\varphi 0}}{\partial A_R} * \Delta^\star v) $
                
        * $  Qjac\_n (V_\varphi,A_R)  = + \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_Z              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0}       * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_R}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0}       * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_R}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0}       * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial A_R}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} * \mu_T * Qvisc\_UZ_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal}             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial A_R}\}_{toroidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial A_R}\}_{toroidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (                + \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal}*\sout{PneoZ}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (\{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal}*V_{Z0}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * (               + \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} * \Delta^\star V_{Z0}) $
                
        * $  Qjac\_k (V_\varphi,A_R)  = + \frac{\partial B_{R0}}{\partial A_R} * \mu_T * Qvisc\_UR_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_R} * \mu_T * Qvisc\_Up_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_\varphi,A_R)  = + \{\frac{\partial B_{Z0}}{\partial A_R}\}_{toroidal} * \mu_T * Qvisc\_UZ_{toroidal, k} $
				
        * $  Qjac\_p (V_\varphi,A_Z)  = + \frac{\partial B_{Z0}}{\partial A_Z} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_Z              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_Z} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_\varphi              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_Z}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_Z}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial A_Z}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{Z0}}{\partial A_Z} * \mu_T * Qvisc\_UZ_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_Z} * \mu_T * Qvisc\_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{poloidal}          $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial A_Z}\}_{poloidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial A_Z}\}_{poloidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (\frac{\partial B_{\varphi 0}}{\partial A_Z}*V_{\varphi 0}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * (              + \frac{\partial B_{Z0}}{\partial A_Z} * \Delta^\star V_{Z0} + \frac{\partial B_{\varphi 0}}{\partial A_Z} * \Delta^\star v) $
                
        * $  Qjac\_n (V_\varphi,A_Z)  = + \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_R              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0}       * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_Z}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0}       * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_Z}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0}       * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial A_Z}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} * \mu_T * Qvisc\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal}             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial A_Z}\}_{toroidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial A_Z}\}_{toroidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (\{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal}*PneoR                  ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (\{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal}*V_{R0}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * (\{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} * \Delta^\star V_{R0} ) $

        * $  Qjac\_k (V_\varphi,A_Z)  = + \frac{\partial B_{Z0}}{\partial A_Z} * \mu_T * Qvisc\_UZ_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_Z} * \mu_T * Qvisc\_Up_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal, k} $
                
        * $  Qjac\_kn(V_\varphi,A_Z)  = + \{\frac{\partial(B_{R0})}{\partial A_Z}\}_{toroidal} * \mu_T * Qvisc\_UR_{toroidal, k} $
				
        * $  Qjac\_p (V_\varphi,A_3)  = + \frac{\partial B_{R 0}}{\partial A_3} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_R              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{Z 0}}{\partial A_3} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_Z              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_3} * \{-v\rho_0\vec{V}_0\cdot \nabla \vec{V}_0-v\rho_0\vec{V}_{dia,0}\cdot\nabla\vec{V}_0\}_\varphi              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_3}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_3}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial A_3}\}_{poloidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{R 0}}{\partial A_3} * \mu_T * Qvisc\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{Z 0}}{\partial A_3} * \mu_T * Qvisc\_UZ_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_3} * \mu_T * Qvisc\_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{poloidal}          $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial A_3}\}_{poloidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial A_3}\}_{poloidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (\frac{\partial B_{R 0}}{\partial A_3}*PneoR    + \frac{\partial B_{Z 0}}{\partial A_3}*\sout{PneoZ}   ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (\frac{\partial B_{R 0}}{\partial A_3}*V_{R0} + \frac{\partial B_{Z 0}}{\partial A_3}*V_{Z0}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * (\frac{\partial B_{R 0}}{\partial A_3} * \Delta^\star V_{R0} + \frac{\partial B_{Z 0}}{\partial A_3} * \Delta^\star V_{Z0} + \frac{\partial B_{\varphi 0}}{\partial A_3} * \Delta^\star v) $
                
        * $  Qjac\_n (V_\varphi,A_3)  = + B_{R0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial A_3}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial A_3}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0}    * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial A_3}\}_{toroidal}        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial A_3}\}_{toroidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial A_3}\}_{toroidal}) $
                
        * $  Qjac\_k (V_\varphi,A_3)  = + \frac{\partial B_{R 0}}{\partial A_3} * \mu_T * Qvisc\_UR_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{Z 0}}{\partial A_3} * \mu_T * Qvisc\_UZ_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + \frac{\partial B_{\varphi 0}}{\partial A_3} * \mu_T * Qvisc\_Up_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + p_0 * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{toroidal, k} $
				
        * $  Qjac\_p (V_\varphi,V_R)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_R}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_R}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial V_R}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0} * \mu_T * Qvisc\_UR\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial V_R}\} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial V_R}\}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (B_{R0}*V_R) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * B_{R0} * \Delta^\star b_f $
                
        * $  Qjac\_n (V_\varphi,V_R)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_R}\}_{toroidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_R}\}_{toroidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial V_R}\}_{toroidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0} * \mu_T * Qvisc\_UR\_UR_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UR_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UR_{toroidal} $
               
        * $  Qjac\_k (V_\varphi,V_R)  = + B_{R0} * \mu_T * Qvisc\_UR\_UR_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UR_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UR_{toroidal, k} $
                
        * $  Qjac\_kn(V_\varphi,V_R)  = + B_{R0} * \mu_T * Qvisc\_UR\_UR_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UR_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UR_{toroidal, kn} $
				
        * $  Qjac\_p (V_\varphi,V_Z)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_Z}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_Z}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial V_Z}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0} * \mu_T * Qvisc\_UR\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UZ_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UZ_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*$\{\frac{\partial \sout{PneoR}}{\partial V_Z}\}$ + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial V_Z}\}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (B_{Z0}*V_Z) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * B_{Z0} * \Delta^\star b_f $

        * $  Qjac\_n (V_\varphi,V_Z)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_Z}\}_{toroidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_Z}\}_{toroidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial V_Z}\}_{toroidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0} * \mu_T * Qvisc\_UR\_UR_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UZ_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UZ_{toroidal} $

        * $  Qjac\_k (V_\varphi,V_Z)  = + B_{R0} * \mu_T * Qvisc\_UR\_UR_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UZ_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UZ_{toroidal, k} $

        * $  Qjac\_kn(V_\varphi,V_Z)  = + B_{R0} * \mu_T * Qvisc\_UR\_UR_{toroidal, kn} $\
                $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_UZ_{toroidal, kn} $\
                $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_UZ_{toroidal, kn} $
				
        * $  Qjac\_p (V_\varphi,V_\varphi)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_\varphi}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_\varphi}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial V_\varphi}\}_{poloidal}           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0} * \mu_T * Qvisc\_UR\_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * S_p * (B_{\varphi 0}*V_\varphi) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \mu_{num} * \Delta^\star v * B_{\varphi 0} * \Delta^\star b_f $

        * $  Qjac\_n (V_\varphi,V_\varphi)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial V_\varphi}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial V_\varphi}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial V_\varphi}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{R0} * \mu_T * Qvisc\_UR\_Up_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_Up_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_Up_{toroidal} $

        * $  Qjac\_k (V_\varphi,V_\varphi)  = + B_{R0} * \mu_T * Qvisc\_UR\_Up_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_Up_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_Up_{toroidal, k} $

        * $  Qjac\_kn(V_\varphi,V_\varphi)  = + B_{R0} * \mu_T * Qvisc\_UR\_Up_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \mu_T * Qvisc\_UZ\_Up_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \mu_T * Qvisc\_Up_Up_{toroidal, kn} $
				
        * $  Qjac\_p (V_\varphi,\rho) = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial \rho}\}_{poloidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial \rho}\}_{poloidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial \rho}\}_{poloidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\rho*(T_{i0}+T_{e0})) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial \rho}\}_{poloidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial \rho}\}_{poloidal}) $

        * $  Qjac\_n (V_\varphi,\rho) = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial \rho}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial \rho}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial \rho}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial \rho}\}_{toroidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial \rho}\}_{toroidal}) $

        * $  Qjac\_k (V_\varphi,\rho) = + (\rho*(T_{i0}+T_{e0})) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} $
				
        * $  Qjac\_p (V_\varphi,T_i)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial T_i}\}_{poloidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial T_i}\}_{poloidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial T_i}\}_{poloidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\rho_0*Ti) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}         $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial T_i}\}_{poloidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial T_i}\}_{poloidal}) $

        * $  Qjac\_n (V_\varphi,T_i)  = + B_{R0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_R}{\partial T_i}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_Z}{\partial T_i}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * -\{\frac{\partial \{\rho_0(\vec{V}_0+\vec{V}_{dia,0})\cdot\nabla \vec{V}_0\}_\varphi}{\partial T_i}\}_{toroidal}              $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - v * (B_{R0}*\{\frac{\partial \sout{PneoR}}{\partial T_i}\}_{toroidal} + B_{Z0}*\{\frac{\partial \sout{PneoZ}}{\partial T_i}\}_{toroidal}) $

        * $  Qjac\_k (V_\varphi,T_i)  = + (\rho_0*Ti) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} $
				
        * $  Qjac\_p (V_\varphi,T_e)  = + B_{R0} * \partial_T \mu * T_e * Qvisc\_UR_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \partial_T \mu * T_e * Qvisc\_UZ_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \partial_T \mu * T_e * Qvisc\_Up_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\rho_0*T_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} $

        * $  Qjac\_k (V_\varphi,T_e)  = + B_{R0} * \partial_T \mu * T_e * Qvisc\_UR_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{Z0} * \partial_T \mu * T_e * Qvisc\_UZ_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + B_{\varphi 0} * \partial_T \mu * T_e * Qvisc\_Up_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\rho_0*T_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} $

    * equation 7   (Density equation)

        * $  Pjac   (\rho,\rho) =   v * \rho $
				
        * $  Qjac\_p (\rho,A_R)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_R}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{poloidal} * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_R}\}_{poloidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $

        * $  Qjac\_n (\rho,A_R)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_R}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal} * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_R}\}_{toroidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $

        * $  Qjac\_k (\rho,A_R)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_R}\}_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal, k} * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_R}\}_{poloidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $

        * $  Qjac\_kn(\rho,A_R)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_R}\}_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_R}\}_{toroidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $
				
        * $  Qjac\_p (\rho,A_Z)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{poloidal} * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_Z}\}_{poloidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $

        * $  Qjac\_n (\rho,A_Z)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal} * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_Z}\}_{toroidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $

        * $  Qjac\_k (\rho,A_Z)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_Z}\}_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal, k} * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_Z}\}_{poloidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $

        * $  Qjac\_kn(\rho,A_Z)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_Z}\}_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_Z}\}_{toroidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * (\vec{B}_0\cdot\nabla \rho_0)       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $
				
        * $  Qjac\_p (\rho,A_3)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_3}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{poloidal} * (\vec{B}_0\cdot\nabla \rho_0)    / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_3}\} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * (\vec{B}_0\cdot\nabla \rho_0)    / B_0^4 * \frac{\partial B_0^2}{\partial A_3} $

        * $  Qjac\_n (\rho,A_3)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_3}\}_{toroidal} $

        * $  Qjac\_k (\rho,A_3)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_3}\}_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{toroidal, k} * (\vec{B}_0\cdot\nabla \rho_0)    / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial A_3}\} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * (\vec{B}_0\cdot\nabla \rho_0)    / B_0^4 * \frac{\partial B_0^2}{\partial A_3} $

        * $  Qjac\_kn(\rho,A_3)  = + \rho_0 * \{\frac{\partial (\vec{V}_{dia,0}\cdot\nabla v)}{\partial A_3}\}_{toroidal, kn} $
				
        * $  Qjac\_p (\rho,V_R)  = - v * ( \rho_0 * \{\nabla\cdot\vec{V}_0\}_{VR} + \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_R} ) $

        * $  Qjac\_p (\rho,V_Z)  = - v * ( \rho_0 * \{\nabla\cdot\vec{V}_0\}_{VZ} + \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_Z} ) $
				
        * $  Qjac\_p (\rho,V_\varphi)  = - v * (                + \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_\varphi} ) $

        * $  Qjac\_n (\rho,V_\varphi)  = - v * ( \rho_0 * \{\nabla\cdot\vec{V}_0\}_{V\varphi, toroidal} ) $
				
        * $  Qjac\_p (\rho,\rho) = - v * ( \rho * divU + \{\frac{\partial\{\vec{V}_0\cdot\nabla \rho\}}{\partial \rho}\}_{poloidal} )                   $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + \rho  * \{\vec{V}_{dia,0}\cdot\nabla v\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial \rho}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - D * \{\nabla\rho\cdot\nabla v\}_{poloidal}                      $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial \rho}\}_{poloidal} / B_0^2 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - D_{\perp, num} * \Delta^\star v * \Delta^\star b_f $

        * $  Qjac\_n (\rho,\rho) = - v * (\{\frac{\partial\{\vec{V}_0\cdot\nabla \rho\}}{\partial \rho}\}_{toroidal} )                                 $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial \rho}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial \rho}\}_{toroidal} / B_0^2 $

        * $  Qjac\_k (\rho,\rho) = + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial \rho}\}_{toroidal, k} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial \rho}\}_{poloidal} / B_0^2 $

        * $  Qjac\_kn(\rho,\rho) = + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial \rho}\}_{toroidal, kn} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - D * \{\nabla\rho\cdot\nabla v\}_{toroidal, kn}                     $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   - (D_\parallel-D) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \{\frac{\partial (\vec{B}_0\cdot\nabla \rho_0)}{\partial \rho}\}_{toroidal} / B_0^2 $
				
        * $  Qjac\_p (\rho,T_i ) = + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial T_i}\}_{poloidal} $

        * $  Qjac\_n (\rho,T_i ) = + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial T_i}\}_{toroidal} $

        * $  Qjac\_k (\rho,T_i ) = + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial T_i}\}_{toroidal, k} $

        * $  Qjac\_kn(\rho,T_i ) = + \rho_0 * \{\frac{\partial{\{\vec{V}_{dia,0}\cdot\nabla v\}}}{\partial T_i}\}_{toroidal, kn} $

    * equation 8   (Ion Temperature  equation)

        * $  Pjac   (T_i,T_i)  =   v * \rho_0 * Ti $

        * $  Pjac   (T_i,\rho) =   v * \rho       * T_{i0} $
				
        * $  Qjac\_p (T_i,A_R)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{poloidal} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_R}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $

        * $  Qjac\_n (T_i,A_R)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_R}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $

        * $  Qjac\_k (T_i,A_R)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_R}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $

        * $  Qjac\_kn(T_i,A_R)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_R}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $
				
        * $  Qjac\_p (T_i,A_Z)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{poloidal} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_Z}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $

        * $  Qjac\_n (T_i,A_Z)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_Z}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $

        * $  Qjac\_k (T_i,A_Z)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_Z}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $

        * $  Qjac\_kn(T_i,A_Z)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_Z}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{i0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $
				
        * $  Qjac\_p (T_i,A_3)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{poloidal} * \vec{B}_0\cdot\nabla T_{i0}    / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_3}\} / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{i0}    / B_0^4 * \frac{\partial B_0^2}{\partial A_3} $

        * $  Qjac\_k (T_i,A_3)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{i0}    / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial A_3}\} / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{i0}    / B_0^4 * \frac{\partial B_0^2}{\partial A_3} $
				
        * $  Qjac\_p (T_i,V_R)  = + v * ( - \rho_0 * \{\vec{V}\cdot\nabla T_i\}_{V_R} - T_{i0} * \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_R}  -  \gamma * p_{i0} * \{\nabla\cdot\vec{V}_0\}_{VR}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_UR_{poloidal} $

        * $  Qjac\_n (T_i,V_R)  = + v * (\gamma-1) * Qvisc\_T\_UR_{toroidal} $
				
        * $  Qjac\_p (T_i,V_Z)  = + v * ( - \rho_0 * \{\vec{V}\cdot\nabla T_i\}_{V_Z} - T_{i0} * \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_Z}  -  \gamma * p_{i0} * \{\nabla\cdot\vec{V}_0\}_{VZ}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_UZ_{poloidal} $

        * $  Qjac\_n (T_i,V_Z)  = + v * (\gamma-1) * Qvisc\_T\_UZ_{toroidal} $
				
        * $  Qjac\_p (T_i,V_\varphi)  = + v * ( - \rho_0 * \{\vec{V}\cdot\nabla T_i\}_{V_\varphi} - T_{i0} * \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_\varphi}                             ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_Up_{poloidal} $

        * $  Qjac\_n (T_i,V_\varphi)  = + v * (                                          -  \gamma * p_{i0} * \{\nabla\cdot\vec{V}_0\}_{V\varphi, toroidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_Up_{toroidal} $
				
        * $  Qjac\_p (T_i,\rho) = + v * ( - \rho * \vec{V}\cdot\nabla T_i - T_{i0} * \{\frac{\partial\{\vec{V}_0\cdot\nabla \rho\}}{\partial \rho}\}_{poloidal}  -  \gamma * (\rho*T_{i0}) * divU) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * \sout{ddTi\_e\_drho} $

        * $  Qjac\_n (T_i,\rho) = + v * (                 - T_{i0} * \{\frac{\partial\{\vec{V}_0\cdot\nabla \rho\}}{\partial \rho}\}_{toroidal}                             ) $
				
        * $  Qjac\_p (T_i,T_i)  = + v * ( - \rho_0 * \frac{\partial (\vec{V}\cdot\nabla T_i)}{\partial T_i}_{poloidal} - Ti * \vec{V}_0\cdot\nabla\rho_0 - \gamma * (\rho_0*Ti) * divU ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_T_{poloidal}                                           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \kappa_i * \{\frac{\partial (\nabla T_i \cdot \nabla v)}{\partial T_i}\}{poloidal}                                         $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial T_i}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\partial_T \kappa_{i,\parallel}*Ti    ) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * \sout{ddTi\_e\_dTi}                                                            $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \kappa_{i,num} * \Delta^\star v * \Delta^\star b_f $

        * $  Qjac\_n (T_i,T_i)  = + v * ( - \rho_0 * \frac{\partial (\vec{V}\cdot\nabla T_i)}{\partial T_i}_{toroidal}                                            ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_T_{toroidal}                                           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial T_i}\}_{toroidal} / B_0^2 $

        * $  Qjac\_k (T_i,T_i)  = - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial T_i}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\partial_T \kappa_{i,\parallel}*Ti    ) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{i0}       / B_0^2 $

        * $  Qjac\_kn(T_i,T_i)  = - \kappa_i * \{\frac{\partial (\nabla T_i \cdot \nabla v)}{\partial T_i}\}{toroidal, kn}                                        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{i,\parallel, T}-\kappa_i) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \{\frac{\partial (\vec{B}\cdot\nabla T_i)}{\partial T_i}\}_{toroidal} / B_0^2 $
				
        * $  Qjac\_p (T_i,T_e)  = + v * \sout{ddTi\_e\_dTe} $

    * equation 9   (Electron Temperature  equation)

        * $  Pjac   (T_e,T_e)  =   v * \rho_0 * T_e $

        * $  Pjac   (T_e,\rho) =   v * \rho       * T_{e0} $
				
        * $  Qjac\_p (T_e,A_R)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{poloidal} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_R}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_R}_{poloidal} $

        * $  Qjac\_n (T_e,A_R)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_R}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_R}_{toroidal} $

        * $  Qjac\_k (T_e,A_R)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_R}\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_R}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{poloidal} $

        * $  Qjac\_kn(T_e,A_R)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_R}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_R}\}_{toroidal} $

        * $  Qjac_pnn(T_e,A_R) = v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_R}_{toroidal,toroidal\ ?,\ what\ is\ nn\ ?} $

        * $  Qjac\_p (T_e,A_Z)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{poloidal} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_Z}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_Z}_{poloidal} $

        * $  Qjac\_n (T_e,A_Z)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_Z}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_Z}_{toroidal} $

        * $  Qjac\_k (T_e,A_Z)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_Z}\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_Z}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{poloidal} $

        * $  Qjac\_kn(T_e,A_Z)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_Z}\}_{toroidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{e0}       / B_0^4 * \{\frac{\partial B_0^2}{\partial A_Z}\}_{toroidal} $

        * $  Qjac_pnn(T_e,A_Z) = v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_Z}_{toroidal,toroidal\ ?,\ what\ is\ nn\ ?} $
				
        * $  Qjac\_p (T_e,A_3)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{poloidal} * \vec{B}_0\cdot\nabla T_{e0}    / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_3}\} / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal}    * \vec{B}_0\cdot\nabla T_{e0}    / B_0^4 * \frac{\partial B_0^2}{\partial A_3} $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_3}_{poloidal} $

        * $  Qjac\_k (T_e,A_3)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\frac{\partial (\vec{B}_0\cdot\nabla v)}{\partial A_3}\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{e0}    / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial A_3}\} / B_0^2             $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k}    * \vec{B}_0\cdot\nabla T_{e0}    / B_0^4 * \frac{\partial B_0^2}{\partial A_3} $

        * $  Qjac\_n (T_e,A_3)  = v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_3}_{toroidal} $

        * $  Qjac_pnn(T_e,A_3) = v * (\gamma-1) * \eta_{T,ohm} * \frac{\partial J^2}{\partial A_3}_{toroidal,toroidal\ ?,\ what\ is\ nn\ ?} $
				
        * $  Qjac\_p (T_e,V_R)  = + v * ( - \rho_0 * \{\vec{V}\cdot\nabla T_e\}_{V_R} - T_{e0} * \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_R} - \gamma * pe0 * \{\nabla\cdot\vec{V}_0\}_{VR}  ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_UR_{poloidal} $

        * $  Qjac\_n (T_e,V_R)  = + v * (\gamma-1) * Qvisc\_T\_UR_{toroidal} $
				
        * $  Qjac\_p (T_e,V_Z)  = + v * ( - \rho_0 * \{\vec{V}\cdot\nabla T_e\}_{V_Z} - T_{e0} * \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_Z} - \gamma * pe0 * \{\nabla\cdot\vec{V}_0\}_{VZ}  ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_UZ_{poloidal} $

        * $  Qjac\_n (T_e,V_Z)  = + v * (\gamma-1) * Qvisc\_T\_UZ_{toroidal} $
				
        * $  Qjac\_p (T_e,V_\varphi)  = + v * ( - \rho_0 * \{\vec{V}\cdot\nabla T_e\}_{V_\varphi} - T_{e0} * \frac{\partial\{\vec{V}\cdot\nabla \rho_0\}}{\partial V_\varphi}                             ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_Up_{poloidal} $

        * $  Qjac\_n (T_e,V_\varphi)  = + v * (                                          -  \gamma * pe0 * \{\nabla\cdot\vec{V}_0\}_{V\varphi, toroidal}) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_Up_{toroidal} $
				
        * $  Qjac\_p (T_e,\rho) = + v * ( - \rho * \vec{V}\cdot\nabla T_e - T_{e0} * \{\frac{\partial\{\vec{V}_0\cdot\nabla \rho\}}{\partial \rho}\}_{poloidal}  - \gamma * (\rho*T_{e0}) * divU ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * \sout{ddTe\_i\_drho} $

        * $  Qjac\_n (T_e,\rho) = + v * (                 - T_{e0} * \{\frac{\partial\{\vec{V}_0\cdot\nabla \rho\}}{\partial \rho}\}_{toroidal}                             ) $
				
        * $  Qjac\_p (T_e,T_e)  = + v * ( - \rho_0 * \frac{\partial (\vec{V}\cdot\nabla T_e)}{\partial T_e}_{poloidal} - T_e * \vec{V}_0\cdot\nabla\rho_0 - \gamma * (\rho_0*T_e) * divU ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_T_{poloidal}                                           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \kappa_e * \{\frac{\partial (\nabla T_e \cdot \nabla v)}{\partial T_e}\}{poloidal}                                         $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial T_e}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\partial_T \kappa_{e,\parallel}*T_e    ) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * \sout{ddTe\_i\_dTe}                                                            $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * \partial_T\eta_{ohm} * T_e * J^2                                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - \kappa_{e,num} * \Delta^\star v * \Delta^\star b_f $

        * $  Qjac\_n (T_e,T_e)  = + v * ( - \rho_0 * \frac{\partial (\vec{V}\cdot\nabla T_e)}{\partial T_e}_{toroidal}                                            ) $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  + v * (\gamma-1) * Qvisc\_T\_T_{toroidal}                                           $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{poloidal} * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial T_e}\}_{toroidal} / B_0^2 $

        * $  Qjac\_k (T_e,T_e)  = - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial T_e}\}_{poloidal} / B_0^2                $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  - (\partial_T \kappa_{e,\parallel}*T_e    ) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \vec{B}_0\cdot\nabla T_{e0}       / B_0^2 $

        * $  Qjac\_kn(T_e,T_e)  = - \kappa_e * \{\frac{\partial (\nabla T_e \cdot \nabla v)}{\partial T_e}\}{toroidal, kn}                                        $\
$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ - (\kappa_{e,\parallel, T}-\kappa_e) * \{\vec{B}_0\cdot\nabla v\}_{toroidal, k} * \{\frac{\partial (\vec{B}\cdot\nabla T_e)}{\partial T_e}\}_{toroidal} / B_0^2 $
				
        * $  Qjac\_p (T_e,T_i)  = + v * \sout{ddTe\_i\_dTi} $
				

    ```fortran
                    enddo ! in loop (n_tor, or not...)
                    enddo ! l loop n_degrees
                    enddo ! k loop (n_vertex)
                enddo ! im loop (n_tor, or not...)
            enddo ! mp loop (n_plane)
            
        enddo ! ms loop
        enddo ! mt loop
    enddo
    enddo
    ```


end subroutine element_matrix_fft
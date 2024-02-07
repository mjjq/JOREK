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
                    PneoR = 0.d0
                    PneoZ = 0.d0
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
    * rhs_p_ij  = 0.d0
    * rhs_k_ij  = 0.d0
    * Pvec_prev = 0.d0 ! The time derivative part
    * Qvec_p    = 0.d0 ! The rest of the RHS (poloidal part)
    * Qvec_k    = 0.d0 ! The rest of the RHS (toroidal part that has phi-derivatives of the test-function)
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
        * $Qvec\_k(\rho)=+\rho_0\vec{V}_0\cdot\nabla v$\
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
                enddo ! im loop (n_tor, or not...)
            enddo ! mp loop (n_plane)
            
        enddo ! ms loop
        enddo ! mt loop
    enddo
    enddo
    ```


end subroutine element_matrix_fft
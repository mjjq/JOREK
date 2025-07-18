#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

e_charge = 1.602176634e-19 # C
mu0      = 4 * np.pi * 1e-7 # Vs/Am

# /u/hazh/SCRATCH/pitagora/20230825_fix_axis_bug_get_all_commits/AUG_36663_3.75_ed22_711_00_ultimate1_eta_zk_fg_q0_gt_1_nmax_4_nflux_80_lft_q0_v4_Dpar0.1_current_source_2_thermalization2_correct_zkperp_etaARAZ_simple
Ti_file = 'jorek_36663_3.75-3.95_Ti_ed22.txt'
Te_file = 'jorek_36663_3.75-3.95_cde_te_ed22.txt'
ne_file = 'jorek_36663_3.75-3.95_ne_from_PT_ed22.txt'
FFp_file = 'jorek_36663_3.75-3.95_FFP_ed22.txt'
central_density = 1.01 * 1e20 # m^-3
psi_axis = 2.092050571E-01
psi_bnd = 2.109960708E-09
R_axis = 1.713940447E+00

# /u/hazh/SCRATCH/pitagora/2025_0502_scan_visco_eta_Hartmann_Prandtl_PhaseII_lower_beta/Base_case_same_q_prof_as_FP_inc_ne_BC_no_particlesource_pitagora5
# Ti_file = 'jorek_36663_2.20-2.40_Ti_ed21.txt'
# Te_file = 'jorek_36663_2.20-2.40_cde_te_ed21.txt'
# ne_file = 'jorek_36663_2.20-2.40_ne_from_PT_ed21.txt'
# FFp_file = 'jorek_36663_2.20-2.40_FFP_ed21_q0_gt_1.txt'
# central_density = 0.905 * 1e20 # m^-3
# psi_axis = 2.230539232E-01
# psi_bnd = -5.124679080E-12
# R_axis = 1.693197427E+00

T_norm = e_charge * mu0 * central_density

data = np.loadtxt(Ti_file)
psin = data[:, 0]
rho_p = psin**0.5
Ti_eV = data[:, 1] / T_norm

data = np.loadtxt(Te_file)
Te_eV = data[:, 1] / T_norm

data = np.loadtxt(ne_file)
ne = data[:,1] * central_density

data = np.loadtxt(FFp_file)
psin_FFp = data[:, 0]
rho_p_FFp = psin_FFp**0.5
FF_prime_FFp = data[:,1]

fq_target = interp1d(psin_FFp, FF_prime_FFp, kind='quadratic', fill_value="extrapolate")
FF_prime = fq_target(psin)

pres = e_charge * ne * (Ti_eV + Te_eV)
psi = psin * (psi_bnd - psi_axis) + psi_axis

p_prime = np.gradient(pres, psi) * mu0 * R_axis **2

plt.plot(rho_p, p_prime, label = 'p_prime')
plt.plot(rho_p, FF_prime, label = 'FF_prime')
plt.plot(rho_p, p_prime + FF_prime, label = 'p_prime + FF_prime')

modify_FFprime_based_on_Pprime = False

if modify_FFprime_based_on_Pprime:

    multiplier = 1.0/3.0
    delta_p_prime = p_prime * (1 - multiplier)

    p_prime_new = p_prime - delta_p_prime
    FF_prime_new = FF_prime - delta_p_prime * 1.0

    data_save=np.column_stack((psin, FF_prime_new))
    file_save = FFp_file[:-4] + '_' + '{:.3f}'.format(multiplier) + '.txt'

    np.savetxt(file_save, data_save, delimiter=' ', comments='')
    print(f"File is saved to {file_save}")

    plt.plot(rho_p, p_prime_new, '--', label = 'p_prime_new')
    plt.plot(rho_p, FF_prime_new, '--', label = 'FF_prime_new')
    plt.plot(rho_p, p_prime_new + FF_prime_new, '--', label = 'p_prime_new + FF_prime_new')

# plt.plot(rho_p, np.gradient(pres, psin), label = 'dpres_dpsin')
plt.legend()
plt.grid()
plt.show()
#! /usr/bin/env python3

#
# Purpose: Truncate the profile of perpendicular heat conductivities in the 
# core region of 1/1 quasi-interchange mode to vary the pressure profile, 
# based on TGLF runs and suggestions from E. Fable and F. Stefanelli.
#
# Date: 2025-05-27
# Author: Haowei Zhang, IPP Garching
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import sys
import math
import plot_multiple_files as pltm

from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.fftpack import fft, ifft

plot_electron = False
if plot_electron:
    file_old_zk_perp = 'jorek_36663_3.75-3.95_kap_e_ed22.txt'
else:
    file_old_zk_perp = 'jorek_36663_3.75-3.95_kap_i_ed22.txt'

file_density     = 'jorek_36663_3.75-3.95_ne_from_PT_ed22.txt'

center_censity   = 1.01e20 # m^-3
mass_proton      = 1.673e-27 # kg
rho0             = center_censity * mass_proton # kg/m3
mu0              = 4 * math.pi * 1e-7 # Vs/Am
gamma            = 5.0/3.0
kap2si           = np.sqrt(rho0 / mu0) / (gamma - 1)

r_truncate       = 0.4
psin_truncate    = r_truncate ** 2
delta_psi        = 0.03
scale_ratio      = 0.1

if __name__ == "__main__":

    column_names, data_zk_perp_old = pltm.read_ascii_file(file_old_zk_perp, 0)
    psin        = data_zk_perp_old[:,0]
    rho_p       = np.sqrt(psin)

    zk_perp_old = data_zk_perp_old[:,1] * kap2si # kg/ms

    column_names_tmp, data_density = pltm.read_ascii_file(file_density, 0)
    data_density[np.where(data_density[:,1] < 5e-2)] = 5e-2
    density =  data_density[:,1] * rho0

    kap_perp_old = zk_perp_old / density # m^2/s
    weight_func = (np.tanh((psin - psin_truncate)/delta_psi) + 1.0) / 2

    # plt.plot(rho_p, zk_perp_old, linestyle='-', label = 'old (kg/ms)')
    plt.plot(rho_p, kap_perp_old, linestyle='-', label = 'old (m$^2$/s)')

    for scale_ratio in [2.5, 2.0, 1.5, 1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 0.8, 0.6, 0.4, 0.2, 0.1]:
        zk_perp_new = zk_perp_old * scale_ratio * (1 - weight_func) + zk_perp_old * weight_func
        kap_perp_new = zk_perp_new / density

        # plt.plot(rho_p, zk_perp_new, linestyle='--', label = fr'new (kg/ms), $\times$ {scale_ratio:.2f}')
        plt.plot(rho_p, kap_perp_new, linestyle='--', label = fr'new (m$^2$/s), $\times$ {scale_ratio:.2f}')

        data_save=np.column_stack((psin, zk_perp_new / kap2si))
        file_save = f"{file_old_zk_perp[:-4]}_x{scale_ratio:.2f}.txt"
        np.savetxt(file_save, data_save, delimiter=' ', comments='')
        print(f"File is saved to {file_save}")
    
    # plt.plot(rho_p, weight_func)
    plt.xlabel('$\sqrt{\psi_n}$')
    if plot_electron:
        plt.ylabel('electron')
    else:
        plt.ylabel('ion')
    plt.grid()
    legend = plt.legend()
    legend.set_draggable(True)
    plt.tight_layout()
    plt.show()

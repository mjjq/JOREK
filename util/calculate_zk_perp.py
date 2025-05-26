#! /usr/bin/env python3

#
# Purpose: Calculate the zk_perp profile from heating source and 
# temperature profile, to satisfy the equilibrium state, 
# that is: div(zk_perp grad T) + S_T = 0
#
# Date: 2025-05-23
# Author: Haowei Zhang, IPP Garching
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import sys
import plot_multiple_files as pltm

from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.fftpack import fft, ifft

def low_pass_filter(data, threshold=0.1):
    f = fft(data)
    f[int(threshold*len(f)):] = 0
    return np.real(ifft(f))

grad_T_flux_file    = 'exprs_surface_flux_s000000.dat'
source_T_int3D_file = 'int3Dprofile_s000000.dat'
gamma = 5.0/3.0
npsi_target = 200

plot_electron = True

if __name__ == "__main__":

    psin_target = np.linspace(0, 1, npsi_target)

    if plot_electron:
        postproc_folder     = 'postproc_electron/'
        file_old_zk_perp = 'jorek_36663_3.75-3.95_kap_e_ed22.txt'
    else:
        postproc_folder     = 'postproc_ion/'
        file_old_zk_perp = 'jorek_36663_3.75-3.95_kap_i_ed22.txt'
    column_names, data_zk_perp_old = pltm.read_ascii_file(file_old_zk_perp, 0)



    column_names, data_grad_T_flux = pltm.read_ascii_file(postproc_folder+grad_T_flux_file, 1)
    if column_names[0] == 'Psi_N':
        psin = data_grad_T_flux[:, 0]
        rho_p = np.sqrt(psin)
    elif column_names[0] == 'rho_p':
        psin = data_grad_T_flux[:, 0] ** 2
        rho_p = data_grad_T_flux[:, 0]
    else:
        sys.exit(f'the data of first column is not clear in {grad_T_flux_file}')
    
    grad_T_flux = data_grad_T_flux[:, 1]
    

    column_names, data_heat_source_int3D = pltm.read_ascii_file(postproc_folder + source_T_int3D_file, 1)
    print(column_names)
    if plot_electron:
        heat_source_int3D = data_heat_source_int3D[:, 4] # electron
    else:
        heat_source_int3D = data_heat_source_int3D[:, 5] # ion

    heat_source_int3D = heat_source_int3D * (gamma - 1)

    zk_perp_original =  - heat_source_int3D / grad_T_flux
    heat_source_int3D = gaussian_filter1d(heat_source_int3D, sigma=1)
    grad_T_flux = gaussian_filter1d(grad_T_flux, sigma=1)
    zk_perp_smoothed =  - heat_source_int3D / grad_T_flux
    zk_perp_smoothed = gaussian_filter1d(zk_perp_smoothed, sigma=1)
    # zk_perp_smoothed = low_pass_filter(zk_perp_smoothed, 4)
    f = interp1d(psin, zk_perp_smoothed, kind='linear', fill_value='extrapolate')
    zk_perp_target = f(psin_target)

    
    plt.plot(rho_p, zk_perp_original, linestyle='-', label = 'unsmoothed')
    plt.plot(rho_p, zk_perp_smoothed, linestyle='--', label = 'smoothed')
    plt.plot(np.sqrt(psin_target), zk_perp_target, linestyle=':', label = 'extrapolated')
    plt.plot(np.sqrt(data_zk_perp_old[:,0]), data_zk_perp_old[:,1], linestyle=':', label = 'old profile')
    # plt.plot(np.sqrt(data_zk_perp_old[:,0]), data_zk_perp_old[:,1] * np.mean(zk_perp_original) / np.mean(data_zk_perp_old[:,1]), linestyle=':', label = 'old profile x ratio')
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

    data_save=np.column_stack((psin_target, zk_perp_target))
    np.savetxt(file_old_zk_perp[:-4]+'_new.txt', data_save, delimiter=' ', comments='')  
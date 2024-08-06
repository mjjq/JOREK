#! /usr/bin/env python3

#
# Purpose: Plot different contributions (dynamo term, resistive term, etc.)
# in induction and current evolution equations.
# Date: 2024-08-06
# Author: Haowei Zhang, IPP Garching
#

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltmp
import plot_poincare_all as pltpoinc
from scipy.ndimage import gaussian_filter1d

folder = 'postproc_v1/'

# balance in base case
pattern0 = os.path.join(folder, 'exprs_averaged_s[3][4-9]???.dat')
pattern1 = os.path.join(folder, 'exprs_averaged_s[4][0-2]???.dat')

# initial stage, dynamo >> resistive, Jphi decreases
# pattern0 = os.path.join(folder, 'exprs_averaged_s[0][0-9]???.dat')
# pattern1 = os.path.join(folder, 'exprs_averaged_s[5][0-2]???.dat')
matched_file_list = [os.path.basename(file) for file in glob.glob(pattern0)] + [os.path.basename(file) for file in glob.glob(pattern1)]
file_list = sorted(matched_file_list, key=lambda x: int(x.split('_s')[-1].split('.dat')[0]))
nfiles = len(file_list)

min_time = pltpoinc.get_time(int(file_list[0].split('_s')[-1].split('.dat')[0]))
max_time = pltpoinc.get_time(int(file_list[-1].split('_s')[-1].split('.dat')[0]))

plot_R2 = False
time2si = 4.607490970E-04

if __name__ == "__main__":

    file_count = 0
    print(f'Total {nfiles} files')

    for file_name in file_list:
        file_path = folder + file_name
        column_names, data = pltmp.read_ascii_file(file_path, skiprows = 1, droprows = -1)
        if file_count == 0:
            print(f'Data list includes {column_names}')
        print(f'Processing No. {file_count} of {nfiles} file: {file_path}')
        
        file_count = file_count + 1

        if file_count == 1:
            psiN = data[:,column_names.index("Psi_N")] / nfiles
            Jtor = data[:,column_names.index("Jtor")] / nfiles

            VxBt_n_eq_one = data[:,column_names.index("VXBt")] / nfiles
            VxBt_n_ge_one = data[:,column_names.index("VXBt_all")] / nfiles
            VxBt_n_eq_zero = data[:,column_names.index("VXBt00")] / nfiles
            etaJtor = data[:,column_names.index("etaJt")] / nfiles
            etaSJtor = data[:,column_names.index("etaSJt")] / nfiles

            if plot_R2:
                VxBt_n_eq_one_R2 = data[:,column_names.index("VXBt_R2")] / nfiles
                VxBt_n_ge_one_R2 = data[:,column_names.index("VXBt_all_R2")] / nfiles
                VxBt_n_eq_zero_R2 = data[:,column_names.index("VXBt00_R2")] / nfiles
                etaJtor_R2 = data[:,column_names.index("etaJt_R2")] / nfiles
                etaSJtor_R2 = data[:,column_names.index("etaSJt_R2")] / nfiles
        else:
            psiN = psiN + data[:,column_names.index("Psi_N")] / nfiles
            Jtor = Jtor + data[:,column_names.index("Jtor")] / nfiles

            VxBt_n_eq_one = VxBt_n_eq_one + data[:,column_names.index("VXBt")] / nfiles
            VxBt_n_ge_one = VxBt_n_ge_one + data[:,column_names.index("VXBt_all")] / nfiles
            VxBt_n_eq_zero = VxBt_n_eq_zero + data[:,column_names.index("VXBt00")] / nfiles
            etaJtor = etaJtor + data[:,column_names.index("etaJt")] / nfiles
            etaSJtor = etaSJtor + data[:,column_names.index("etaSJt")] / nfiles

            if plot_R2:
                VxBt_n_eq_one_R2 = VxBt_n_eq_one_R2 + data[:,column_names.index("VXBt_R2")] / nfiles
                VxBt_n_ge_one_R2 = VxBt_n_ge_one_R2 + data[:,column_names.index("VXBt_all_R2")] / nfiles
                VxBt_n_eq_zero_R2 = VxBt_n_eq_zero_R2 + data[:,column_names.index("VXBt00_R2")] / nfiles
                etaJtor_R2 = etaJtor_R2 + data[:,column_names.index("etaJt_R2")] / nfiles
                etaSJtor_R2 = etaSJtor_R2 + data[:,column_names.index("etaSJt_R2")] / nfiles

    VxBt_n_all = VxBt_n_ge_one + VxBt_n_eq_zero
    VxBt_n_all_origin = VxBt_n_all
    etaJtor_origin = etaJtor
    etaSJtor_origin = etaSJtor

    sigma_VxB = 3
    sigma_etaJ = 3
    VxBt_n_all = gaussian_filter1d(VxBt_n_all, sigma_VxB)
    VxBt_n_eq_one = gaussian_filter1d(VxBt_n_eq_one, sigma_VxB)
    VxBt_n_ge_one = gaussian_filter1d(VxBt_n_ge_one, sigma_VxB)
    VxBt_n_eq_zero = gaussian_filter1d(VxBt_n_eq_zero, sigma_VxB)
    etaJtor = gaussian_filter1d(etaJtor, sigma_etaJ)
    etaSJtor = gaussian_filter1d(etaSJtor, sigma_etaJ)

    rhop = np.sqrt(psiN)
    eta_T = etaJtor / Jtor


    ddrho_VxBt_n_all = np.gradient(VxBt_n_all, rhop)
    d2drho2_VxBt_n_all = np.gradient(ddrho_VxBt_n_all, rhop)

    ddrho_VxBt_n_ge_one = np.gradient(VxBt_n_ge_one, rhop)
    d2drho2_VxBt_n_ge_one = np.gradient(ddrho_VxBt_n_ge_one, rhop)

    ddrho_VxBt_n_eq_zero = np.gradient(VxBt_n_eq_zero, rhop)
    d2drho2_VxBt_n_eq_zero = np.gradient(ddrho_VxBt_n_eq_zero, rhop)

    ddrho_etaJtor = np.gradient(etaJtor, rhop)
    d2drho2_etaJtor = np.gradient(ddrho_etaJtor, rhop)

    ddrho_etaSJtor = np.gradient(etaSJtor, rhop)
    d2drho2_etaSJtor = np.gradient(ddrho_etaSJtor, rhop)

    # testf = np.sin(4*rhop) + np.cos(8*rhop)
    # ddrho_testf_ana = 4*np.cos(4*rhop)-8*np.sin(8*rhop)
    # d2drho2_testf_ana = -16*np.sin(4*rhop)-64*np.cos(8*rhop)
    # ddrho_testf = np.gradient(testf, rhop)
    # d2drho2_testf = np.gradient(ddrho_testf, rhop)

    # plt.figure()
    # plt.plot(rhop, d2drho2_testf_ana, 'o', linestyle = '--', label = 'Analytical solution')
    # plt.plot(rhop, d2drho2_testf, 's', linestyle = ':', label = 'Numerical solution')
    # plt.xlabel(r'$\sqrt{\psi_N}$')
    # plt.legend()
    # plt.grid()


    # plt.figure()
    # plt.semilogy(rhop, eta_T)
    # plt.xlabel(r'$\sqrt{\psi_N}$')
    # plt.ylabel(r'$\eta_{spz}$')
    # plt.grid()

    plt.figure()
    plt.plot(rhop, Jtor, label = r'$J_\varphi$')
    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'Different terms')
    plt.legend()
    plt.grid()
    plt.title(f'Averaged between {"{:.1f}".format(min_time * time2si)} and {"{:.1f}".format(max_time * time2si)} [ms]')

    if plot_R2:
        plt.figure()
        plt.plot(rhop, VxBt_n_ge_one_R2 + VxBt_n_eq_zero_R2, label = r'$(VxB_{n>0} + VxB_{n=0})/R^2$')
        plt.plot(rhop, -etaJtor_R2 + etaSJtor_R2, label = r'$-\eta (J_{\varphi} - S_J)/R^2$')
        plt.xlabel(r'$\sqrt{\psi_N}$')
        plt.ylabel(r'Different terms')
        plt.legend()
        plt.grid()
        plt.title(f'Averaged between {"{:.1f}".format(min_time * time2si)} and {"{:.1f}".format(max_time * time2si)} [ms]')

    plt.figure()
    plt.plot(rhop, VxBt_n_ge_one + VxBt_n_eq_zero, color='red', label = r'$(VxB)_{all}=(VxB)_{n>0} + (VxB)_{n=0}$')
    plt.plot(rhop, -etaJtor + etaSJtor, color='green', label = r'$-\eta (J_{\varphi} - S_J)$')

    plt.plot(rhop, VxBt_n_all_origin, 'r.',linestyle = 'none' ,label = r'$(VxB)_{all, unsmoothed}=(VxB)_{n>0} + (VxB)_{n=0}$')
    plt.plot(rhop, -etaJtor_origin + etaSJtor_origin, 'gx',linestyle = 'none', label = r'$-\eta (J_{\varphi} - S_J)_{unsmoothed}$')

    plt.plot(rhop, VxBt_n_ge_one, '.', linestyle = '--', label = r'$VxB_{n>0}$')
    plt.plot(rhop, VxBt_n_eq_zero, linestyle = 'dotted', label = r'$VxB_{n=0}$')

    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'Different terms')
    plt.legend()
    plt.grid()
    plt.title(f'Averaged between {"{:.1f}".format(min_time * time2si)} and {"{:.1f}".format(max_time * time2si)} [ms]')

    plt.figure()
    # plt.plot(rhop, -d2drho2_VxBt_n_ge_one - d2drho2_VxBt_n_eq_zero, label = r'$-\partial_\rho^2[(VxB)_{n>0} + (VxB)_{n=0}]$')
    # plt.plot(rhop, -d2drho2_VxBt_n_ge_one, label = r'$-\partial_\rho^2 (VxB)_{n>0}$')
    # plt.plot(rhop, -d2drho2_VxBt_n_eq_zero, label = r'$-\partial_\rho^2 (VxB)_{n=0}$')
    plt.plot(rhop, -d2drho2_VxBt_n_all, label = r'$-\partial_\rho^2 (VxB)_{all}$')
    plt.plot(rhop, d2drho2_etaJtor - d2drho2_etaSJtor, label = r'$\partial_\rho^2[\eta (J_{\varphi} - S_J)]$')
    # plt.plot(rhop, VxBt_n_all)
    # plt.plot(rhop, gaussian_filter1d(VxBt_n_all,3), '.')

    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'Different terms')
    legend = plt.legend()  # Show legends for each line
    legend.set_draggable(True)
    plt.grid()
    # plt.tight_layout()
    plt.title(f'Averaged between {"{:.1f}".format(min_time * time2si)} and {"{:.1f}".format(max_time * time2si)} [ms]')

    plt.show()



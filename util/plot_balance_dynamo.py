#! /usr/bin/env python3

#
# Purpose: Plot different contributions (dynamo term, resistive term, etc.)
# in induction and current evolution equations.
# Date: 2024-08-06
# Author: Haowei Zhang, IPP Garching
#

import glob
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltmp
import plot_poincare_all as pltpoinc
from scipy.ndimage import gaussian_filter1d

# Set global font properties using rcParams
plt.rcParams.update({
    'axes.formatter.limits': (-3, 3),
    'axes.formatter.useoffset': True,  # Use scientific notation
    'font.size': 14,                # Global font size
    'lines.linewidth': 2.0,         # Global line width
    'axes.linewidth': 2.0           # Set line width for both x and y axes
    # 'font.family': 'serif',         # Use serif fonts, which includes Times New Roman
    # 'font.serif': ['Times New Roman'],  # Specifically set Times New Roman
    # 'font.weight': 'bold'           # Font weight for all text elements
    # 'legend.fontsize': 14,          # Font size for the legend
    # 'axes.labelsize': 14,           # Font size for axis labels
    # 'axes.titlesize': 16,           # Font size for title
    # 'xtick.labelsize': 12,          # Font size for x-tick labels
    # 'ytick.labelsize': 12,          # Font size for y-tick labels
})

# load the 3D data
folder = 'postproc/'
filenamerule0 = 'exprs_averaged_s00000.dat' # initial
filenamerule1 = 'exprs_averaged_s[3][6-9]???.dat'
filenamerule2 = 'exprs_averaged_s[4][0-2]???.dat'

# initial stage, dynamo >> resistive, Jphi decreases
# filenamerule1 = 'exprs_averaged_s[0][3-5]???.dat'
# filenamerule2 = 'exprs_averaged_s[9][0-2]???.dat'

pattern1 = os.path.join(folder, filenamerule1)
pattern2 = os.path.join(folder, filenamerule2)

file0 = folder + filenamerule0
matched_file_list = [os.path.basename(file) for file in glob.glob(pattern1)] + [os.path.basename(file) for file in glob.glob(pattern2)]
file_list = sorted(matched_file_list, key=lambda x: int(x.split('_s')[-1].split('.dat')[0]))
nfiles = len(file_list)

min_time = pltpoinc.get_time(int(file_list[0].split('_s')[-1].split('.dat')[0]))
max_time = pltpoinc.get_time(int(file_list[-1].split('_s')[-1].split('.dat')[0]))

# load the 2D data
folder2D = '../Base_case2D/postproc/'
filenamerule2D = 'exprs_averaged_s[3][3]????.dat'
pattern2D = os.path.join(folder2D, filenamerule2D)
matched_file_list2D = [os.path.basename(file) for file in glob.glob(pattern2D)]
file_list2D = sorted(matched_file_list2D, key=lambda x: int(x.split('_s')[-1].split('.dat')[0]))
nfiles2D = len(file_list2D)
min_time2D = pltpoinc.get_time(int(file_list2D[0].split('_s')[-1].split('.dat')[0]), '../Base_case2D/')
max_time2D = pltpoinc.get_time(int(file_list2D[-1].split('_s')[-1].split('.dat')[0]), '../Base_case2D/')

plot_R2 = False

time2si         = 4.607490970E-04  # to [ms]
MU_zero         = 4.e-7*math.pi
central_density = 1.01
central_mass    = 1.0
mass_proton     = 1.67262178e-27
rho_norm        = central_density *1.e20 * central_mass * mass_proton
fact_J_mu_zero  = MU_zero
fact_time       = 1/4.607490970E-04
fact_eta        = 1.0 / math.sqrt ( MU_zero / rho_norm )
fact_VXB        = math.sqrt( MU_zero * central_mass * central_density * 1.e20 * mass_proton )

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
            gradPsi2 = data[:,column_names.index("gradPsi2")] / nfiles
            lapPsi = data[:,column_names.index("lapPsi")] / nfiles

            # read the initial equilibrium
            column_names0, data0 = pltmp.read_ascii_file(file0, skiprows = 1, droprows = -1)
            Jtor0 = data0[:,column_names0.index("Jtor")]
            etaJtor0 = data0[:,column_names0.index("etaJt")]
            etaSJtor0 = data0[:,column_names0.index("etaSJt")]
            eta_T0 = etaJtor0 / Jtor0

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
            gradPsi2 = gradPsi2 + data[:,column_names.index("gradPsi2")] / nfiles
            lapPsi = lapPsi + data[:,column_names.index("lapPsi")] / nfiles

            if plot_R2:
                VxBt_n_eq_one_R2 = VxBt_n_eq_one_R2 + data[:,column_names.index("VXBt_R2")] / nfiles
                VxBt_n_ge_one_R2 = VxBt_n_ge_one_R2 + data[:,column_names.index("VXBt_all_R2")] / nfiles
                VxBt_n_eq_zero_R2 = VxBt_n_eq_zero_R2 + data[:,column_names.index("VXBt00_R2")] / nfiles
                etaJtor_R2 = etaJtor_R2 + data[:,column_names.index("etaJt_R2")] / nfiles
                etaSJtor_R2 = etaSJtor_R2 + data[:,column_names.index("etaSJt_R2")] / nfiles

    file_count2D = 0
    print(f'Total {nfiles2D} files')
    for file_name2D in file_list2D:
        file_path2D = folder2D + file_name2D
        column_names2D, data2D = pltmp.read_ascii_file(file_path2D, skiprows = 1, droprows = -1)
        if file_count2D == 0:
            print(f'Data list includes {column_names2D}')
        print(f'Processing No. {file_count2D} of {nfiles2D} file: {file_path2D}')

        file_count2D = file_count2D + 1

        if file_count2D == 1:
            psiN2D = data2D[:,column_names2D.index("Psi_N")] / nfiles2D
            Jtor2D = data2D[:,column_names2D.index("Jtor")] / nfiles2D
            etaJtor2D = data2D[:,column_names2D.index("etaJt")] / nfiles2D
            etaSJtor2D = data2D[:,column_names2D.index("etaSJt")] / nfiles2D

        else:
            psiN2D = psiN2D + data2D[:,column_names2D.index("Psi_N")] / nfiles2D
            Jtor2D = Jtor2D + data2D[:,column_names2D.index("Jtor")] / nfiles2D
            etaJtor2D = etaJtor2D + data2D[:,column_names2D.index("etaJt")] / nfiles2D
            etaSJtor2D = etaSJtor2D + data2D[:,column_names2D.index("etaSJt")] / nfiles2D

        eta_T2D = etaJtor2D / Jtor2D


    VxBt_n_all = VxBt_n_ge_one + VxBt_n_eq_zero
    VxBt_n_all_origin = VxBt_n_all
    etaJtor_origin = etaJtor
    etaSJtor_origin = etaSJtor

    sigma_VxB = 3 #* 2.5 # = 3 for 200 points; = 3 * 2.5 for 500 points
    sigma_etaJ = 3 #* 2.5
    VxBt_n_all = gaussian_filter1d(VxBt_n_all, sigma_VxB)
    VxBt_n_eq_one = gaussian_filter1d(VxBt_n_eq_one, sigma_VxB)
    VxBt_n_ge_one = gaussian_filter1d(VxBt_n_ge_one, sigma_VxB)
    VxBt_n_eq_zero = gaussian_filter1d(VxBt_n_eq_zero, sigma_VxB)
    etaJtor = gaussian_filter1d(etaJtor, sigma_etaJ)
    etaSJtor = gaussian_filter1d(etaSJtor, sigma_etaJ)
    Jtor = gaussian_filter1d(Jtor, sigma_etaJ)

    rhop = np.sqrt(psiN)
    eta_T = etaJtor / Jtor
    SJtor = etaSJtor / eta_T


    ddrho_VxBt_n_all = np.gradient(VxBt_n_all, psiN)
    d2drho2_VxBt_n_all = np.gradient(ddrho_VxBt_n_all, psiN) * gradPsi2 + ddrho_VxBt_n_all * lapPsi

    ddrho_VxBt_n_ge_one = np.gradient(VxBt_n_ge_one, psiN)
    d2drho2_VxBt_n_ge_one = np.gradient(ddrho_VxBt_n_ge_one, psiN) * gradPsi2 + ddrho_VxBt_n_ge_one * lapPsi

    ddrho_VxBt_n_eq_zero = np.gradient(VxBt_n_eq_zero, psiN)
    d2drho2_VxBt_n_eq_zero = np.gradient(ddrho_VxBt_n_eq_zero, psiN) * gradPsi2 + ddrho_VxBt_n_eq_zero * lapPsi

    ddrho_etaJtor = np.gradient(etaJtor, psiN)
    d2drho2_etaJtor = np.gradient(ddrho_etaJtor, psiN) * gradPsi2 + ddrho_etaJtor * lapPsi

    ddrho_etaSJtor = np.gradient(etaSJtor, psiN)
    d2drho2_etaSJtor = np.gradient(ddrho_etaSJtor, psiN) * gradPsi2 + ddrho_etaSJtor * lapPsi

    plt.figure()
    plt.semilogy(rhop, eta_T0 / fact_eta, label = r'$\eta_{spz}^{initial}\ [\Omega m]$')
    plt.semilogy(rhop, eta_T / fact_eta, label = r'$\eta_{spz}\ [\Omega m]$')
    plt.semilogy(rhop, eta_T2D / fact_eta, label = r'$\eta_{spz, 2D}\ [\Omega m]$')
    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'$\eta_{spz}$')
    legend = plt.legend()
    legend.set_draggable(True)
    plt.grid()
    plt.title(f'Average: 3D {"{:.0f}".format(min_time * time2si)}~{"{:.0f}".format(max_time * time2si)}; 2D {"{:.0f}".format(min_time2D * time2si)}~{"{:.0f}".format(max_time2D * time2si)} [ms]')
    plt.tight_layout()

    plt.figure()
    plt.plot(rhop, Jtor0/fact_J_mu_zero, label = r'$J_\varphi^{initial}\ [Am^{-2}]$')
    plt.plot(rhop, Jtor/fact_J_mu_zero, label = r'$J_\varphi\ [Am^{-2}]$')
    plt.plot(rhop, Jtor2D/fact_J_mu_zero, label = r'$J_{\varphi, 2D}\ [Am^{-2}]$')
    plt.plot(rhop, SJtor/fact_J_mu_zero, label = r'$S_{j,\varphi}\ [Am^{-2}]$')
    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'different terms')
    legend = plt.legend()
    legend.set_draggable(True)
    # plt.gca().ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.grid()
    plt.title(f'Average: 3D {"{:.0f}".format(min_time * time2si)}~{"{:.0f}".format(max_time * time2si)}; 2D {"{:.0f}".format(min_time2D * time2si)}~{"{:.0f}".format(max_time2D * time2si)} [ms]')
    plt.tight_layout()

    if plot_R2:
        plt.figure()
        plt.plot(rhop, VxBt_n_ge_one_R2 + VxBt_n_eq_zero_R2, label = r'$(VxB_{n>0} + VxB_{n=0})/R^2$')
        plt.plot(rhop, -etaJtor_R2 + etaSJtor_R2, label = r'$-\eta (J_{\varphi} - S_J)/R^2$')
        plt.xlabel(r'$\sqrt{\psi_N}$')
        plt.ylabel(r'Different terms')
        legend = plt.legend()
        legend.set_draggable(True)
        plt.grid()
        plt.title(f'Average: 3D {"{:.0f}".format(min_time * time2si)}~{"{:.0f}".format(max_time * time2si)} [ms]')
        plt.tight_layout()

    plt.figure()
    plt.plot(rhop, VxBt_n_eq_zero / fact_VXB, linestyle = 'dotted', color = 'C1', label = r'$V_{n=0}xB_{n=0}\ [Vm^{-1}]$')
    plt.plot(rhop, VxBt_n_ge_one / fact_VXB, linestyle = '--', color = 'C0', label = r'$V_{n>0}xB_{n>0}\ [Vm^{-1}]$')
    
    plt.plot(rhop, (VxBt_n_ge_one + VxBt_n_eq_zero) / fact_VXB, color='red', label = r'$(V_{all}xB_{all})\ [Vm^{-1}]$')
    plt.plot(rhop, (-etaJtor + etaSJtor) / fact_eta / fact_J_mu_zero, color='green', label = r'$-\eta (J_{\varphi} - S_J)\ [Vm^{-1}]$')

    plt.plot(rhop, -eta_T2D * (Jtor - Jtor2D) / fact_eta / fact_J_mu_zero, color='black', label = r'$-\eta_{2D}\Delta J_\varphi\ [Vm^{-1}]$')
    plt.plot(rhop, -(eta_T - eta_T2D) * Jtor2D / fact_eta / fact_J_mu_zero, color='blue', label = r'$-\Delta\eta J_{\varphi, 2D}\ [Vm^{-1}]$')

    plt.plot(rhop, VxBt_n_all_origin / fact_VXB, 'r.',linestyle = 'none')# ,label = r'$(VxB)_{all, unsmoothed}\ [Vm^{-1}]$')
    plt.plot(rhop, (-etaJtor_origin + etaSJtor_origin) / fact_eta / fact_J_mu_zero, 'gx',linestyle = 'none')#, label = r'$-\eta (J_{\varphi} - S_J)_{unsmoothed}\ [Vm^{-1}]$')


    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'Different terms')
    legend = plt.legend()
    legend.set_draggable(True)
    plt.grid()
    plt.title(f'Average: 3D {"{:.0f}".format(min_time * time2si)}~{"{:.0f}".format(max_time * time2si)}; 2D {"{:.0f}".format(min_time2D * time2si)}~{"{:.0f}".format(max_time2D * time2si)} [ms]')
    plt.tight_layout()

    plt.figure()
    # plt.plot(rhop, -d2drho2_VxBt_n_ge_one - d2drho2_VxBt_n_eq_zero, label = r'$-\partial_\rho^2[(VxB)_{n>0} + (VxB)_{n=0}]$')
    # plt.plot(rhop, -d2drho2_VxBt_n_ge_one, label = r'$-\partial_\rho^2 (VxB)_{n>0}$')
    # plt.plot(rhop, -d2drho2_VxBt_n_eq_zero, label = r'$-\partial_\rho^2 (VxB)_{n=0}$')
    print(f'fact_eta = {fact_eta}')
    plt.plot(rhop, -d2drho2_VxBt_n_all / fact_VXB / MU_zero, label = r'$-\langle\nabla^2 (VxB)_{\varphi,all}\rangle/\mu_0\ [Am^{-2}s^{-1}]$')
    plt.plot(rhop, (d2drho2_etaJtor - d2drho2_etaSJtor) / fact_eta / fact_J_mu_zero / MU_zero, label = r'$\langle\nabla^2[\eta (J_{\varphi} - S_{j,\varphi})]\rangle/\mu_0\ [Am^{-2}s^{-1}]$')
    # plt.plot(rhop, VxBt_n_all)
    # plt.plot(rhop, gaussian_filter1d(VxBt_n_all,3), '.')

    plt.xlabel(r'$\sqrt{\psi_N}$')
    plt.ylabel(r'Different terms')
    legend = plt.legend()
    legend.set_draggable(True)
    plt.grid()
    plt.title(f'Average: 3D {"{:.0f}".format(min_time * time2si)}~{"{:.0f}".format(max_time * time2si)} [ms]')
    plt.tight_layout()

    plt.show()


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
    # plt.tight_layout()
    # plt.show()



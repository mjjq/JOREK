#! /usr/bin/env python3

#
# Purpose: Calculate the target F and FFprime with specified q-profile (unfinished).
#
# Date: 2024-01-27
# Author: Haowei Zhang, IPP Garching
#
# import os
# import math
# import argparse
# import copy
# import plot_vtk as pltvtk
# import plot_q_versus_time as pltqvst
import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltmpf
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d


interp_method = 'linear' # 'linear', 'cubic'

if __name__ == "__main__":
    target_q_file  = 'q_profile_target_Isabel.dat'
    current_q_file = 'q_profile_hw.dat'
    
    output_psi_n = np.linspace(0, 1.1, 221)

    column_names, target_data = pltmpf.read_ascii_file(target_q_file, 1)
    target_psi_n = target_data[:,0]
    target_q = target_data[:,1]
    interp_func_q = interp1d(target_psi_n, target_q, kind=interp_method, fill_value='extrapolate')

    interp_target_q = interp_func_q(output_psi_n)

    # plt.plot(np.sqrt(target_psi_n), abs(target_q))
    # plt.plot(np.sqrt(output_psi_n), abs(interp_target_q))
    # plt.show()

    column_names, current_data = pltmpf.read_ascii_file(current_q_file, 1)
    current_psi_n = current_data[:,0]
    current_q = current_data[:,1]
    current_B_p_intg = current_data[:,3]
    current_psi = current_data[:,4]
    current_F = current_data[:,5]

    interp_func_q2 = interp1d(current_psi_n, current_q, kind=interp_method, fill_value='extrapolate')
    interp_current_q = interp_func_q2(output_psi_n)

    interp_func_Bp_intg = interp1d(current_psi_n, current_B_p_intg, kind=interp_method, fill_value='extrapolate')
    interp_Bp_intg = interp_func_Bp_intg(output_psi_n)

    interp_func_F = interp1d(current_psi_n, current_F, kind=interp_method, fill_value='extrapolate')
    interp_F = interp_func_F(output_psi_n)

    interp_func_psi = interp1d(current_psi_n, current_psi, kind=interp_method, fill_value='extrapolate')
    interp_psi = interp_func_psi(output_psi_n)
    interp_FFp = interp_F * np.gradient(interp_F, interp_psi)

    # plt.plot(np.sqrt(current_psi_n), abs(current_B_p_intg))
    # plt.plot(np.sqrt(output_psi_n), abs(interp_Bp_intg))
    # plt.plot(np.sqrt(output_psi_n), -interp_FFp)
    # plt.grid()
    # plt.show()

    factor=1.0
    interp_F_profile = (interp_target_q * (1.0-factor) + interp_current_q * factor ) / interp_Bp_intg
    interp_F_profile = savgol_filter(interp_F_profile, 21, 1)  # 3 is the order of the filter

    b=0.4
    # b=0.4
    c=0.08
    fq=1-(np.tanh((np.sqrt(output_psi_n)-b)/c)+1)/2.0

    b=0.7
    # b=0.4
    c=0.08
    f=1-(np.tanh((np.sqrt(output_psi_n)-b)/c)+1)/2.0
    delta_q = interp_target_q - interp_current_q

    interp_F_profile = interp_F - (interp_F_profile - interp_F) * f * 0.275 * delta_q/delta_q[0]*fq
    
    # interp_F_profile = interp_F - (interp_F_profile - interp_F) * f * 0.7
    interp_F_profile[0] = interp_F_profile[0] - 0.00001*0.41
    # interp_F_profile = interp_F - (interp_target_q - interp_current_q) * f * 0.01
    plt.plot(np.sqrt(output_psi_n), delta_q*fq/delta_q[0])
    plt.grid()
    plt.show()    
    # interp_F_profile = gaussian_filter1d(interp_F_profile, sigma=5)
    interp_dFdPsi_profile = np.gradient(interp_F_profile, interp_psi)
    interp_FFprime_profile = interp_F_profile * interp_dFdPsi_profile
    # interp_FFprime_profile = gaussian_filter1d(interp_FFprime_profile, sigma=5)
    # plt.plot(np.sqrt(output_psi_n), -interp_current_q)
    # plt.plot(np.sqrt(output_psi_n), -interp_target_q)
    plt.plot(np.sqrt(output_psi_n), interp_F_profile, 'o')
    plt.plot(np.sqrt(output_psi_n), interp_F)
    plt.grid()
    plt.show()    

    interp_FFprime_profile[201:]=0
    # interp_FFprime_profile[3] = interp_FFprime_profile[3] - 0.03
    # interp_FFprime_profile[4] = interp_FFprime_profile[4] - 0.02
    plt.figure()
    plt.plot(np.sqrt(output_psi_n), -interp_FFprime_profile, 'o')
    plt.plot(np.sqrt(output_psi_n), -interp_FFp)
    
    plt.grid()
    plt.show()

    data_save=np.column_stack((output_psi_n, -interp_FFprime_profile))
    np.savetxt('input_profiles/ffprime_jorek.dat', data_save, delimiter=' ', comments='')  
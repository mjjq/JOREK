#! /usr/bin/env python3

# Read files from jorek2_postproc and generate inputs for HELENA

import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltm
from scipy.interpolate import interp1d

numstep = 90000
str_num = f"{numstep:06d}"
postproc_folder = 'postproc/'
postproc_qfile = 'qprofile_s' + str_num + '.dat'
postproc_exprsfile = 'exprs_averaged_s' + str_num + '.dat'
postproc_intg3Dfile = 'int3Dprofile_s' + str_num + '.dat'
postproc_RZfile = 'RZ_maxmin_profile_s' + str_num + '.dat'
e_charge = 1.602176634e-19 #C
# Iptot =  8.129874792E+05 # A
# psi_boundary = 0.00004 # 000300
psi_boundary = -0.00192 # 090000
# psi_axis = 0.20921
# R_axis = 1.713934595E+00
# Z_axis = 7.092228919E-02

psi_n_boundary = [
    1.005025125628141E-01,
    2.010050251256282E-01,
    3.015075376884422E-01,
    4.020100502512563E-01,
    4.974874371859297E-01,
    5.979899497487436E-01,
    6.984924623115578E-01,
    7.989949748743720E-01,
    8.994974874371859E-01,
    9.949748743718594E-01]

# q_boundary = [ # 000300
#     -1.026065640528196E+00,
#     -1.068303009933490E+00,
#     -1.189020649092972E+00,
#     -1.406867805108305E+00,
#     -1.699059254213150E+00,
#     -2.118441230950786E+00,
#     -2.700298699012887E+00,
#     -3.523405034536250E+00,
#     -4.665270925583930E+00,
#     -6.309133434135009E+00]

q_boundary = [  # 090000
    -1.009989030708960E+00,
    -1.068700182437862E+00,
    -1.199657552771730E+00,
    -1.425876174337616E+00,
    -1.721632483593031E+00,
    -2.143848085704047E+00,
    -2.737888838171731E+00,
    -3.488829465250061E+00,
    -4.340120016263127E+00,
    -7.033335088784463E+00]

if __name__ == "__main__":

    column_names, axis_data = pltm.read_ascii_file('axis.dat', 1, droprows = -1) # obtained by plot_live_data -q axis
    R_axis   = axis_data[numstep - 1, 1]
    Z_axis   = axis_data[numstep - 1, 2]
    psi_axis = axis_data[numstep - 1, 3]

    column_names, Ip_data = pltm.read_ascii_file('current.dat', 1, droprows = -1) # obtained by plot_live_data -q cur
    Iptot = axis_data[numstep - 1, 1]

    column_names, data_q = pltm.read_ascii_file(postproc_folder + postproc_qfile, 0, droprows = -1)
    psin_q = data_q[:, 0]
    q = data_q[:, 1]

    column_names_exprs, data_exprs = pltm.read_ascii_file(postproc_folder + postproc_exprsfile, 1, droprows = -2)
    psin_exprs = data_exprs[:, 0]
    ne = data_exprs[:, 1]
    Ti = data_exprs[:, 2]
    Te = data_exprs[:, 3]
    p = data_exprs[:, 4]
    Jtor = data_exprs[:, 8]
    psi = data_exprs[:, 38]

    column_names_intg3D, data_intg3D = pltm.read_ascii_file(postproc_folder + postproc_intg3Dfile, 1, droprows = -1)
    psin_intg3D = data_intg3D[:, 0]
    volume = data_intg3D[:, 1]
    area = data_intg3D[:, 2]
    Ip = data_intg3D[:, 3]

    column_names_RZ, data_RZ = pltm.read_ascii_file(postproc_folder + postproc_RZfile, 1, droprows = -1)
    psin_RZ = data_RZ[:, 0]
    Rmax = data_RZ[:, 1]
    Rmin = data_RZ[:, 2]
    Zmax = data_RZ[:, 3]
    Zmin = data_RZ[:, 4]
    a_minor = (Rmax - Rmin) / 2 # minor radius

    n_target   = 200
    psin_target = np.linspace(0, 1, n_target)
    # q_target = np.interp(psin_target, psin_q, q)
    fq_target = interp1d(psin_q, q, kind='quadratic', fill_value="extrapolate")
    q_target = fq_target(psin_target)

    fne_target = interp1d(psin_exprs, ne, kind='quadratic', fill_value="extrapolate")
    ne_target = fne_target(psin_target)

    fTi_target = interp1d(psin_exprs, Ti, kind='quadratic', fill_value="extrapolate")
    Ti_target = fTi_target(psin_target)

    fTe_target = interp1d(psin_exprs, Te, kind='quadratic', fill_value="extrapolate")
    Te_target = fTe_target(psin_target)

    # fp_target = interp1d(psin_exprs, p, kind='linear', fill_value="extrapolate")
    # p_target = fp_target(psin_target)
    p_target = e_charge * ne_target * (Ti_target + Te_target)
    # if p_target[-1] < 0:
    #     p_target[-1] = 0


    fJtor_target = interp1d(psin_exprs, Jtor, kind='linear', fill_value="extrapolate")
    Jtor_target = fJtor_target(psin_target)

    fpsi_target = interp1d(psin_exprs, psi, kind='linear', fill_value="extrapolate")
    psi_target = fpsi_target(psin_target)

    fvolume_target = interp1d(psin_intg3D, volume, kind='linear', fill_value="extrapolate")
    volume_target = fvolume_target(psin_target)

    farea_target = interp1d(psin_intg3D, area, kind='linear', fill_value="extrapolate")
    area_target = farea_target(psin_target)

    fIp_target = interp1d(psin_intg3D, Ip, kind='linear', fill_value="extrapolate")
    Ip_target = fIp_target(psin_target)

    fa_target = interp1d(psin_RZ, a_minor, kind='linear', fill_value="extrapolate")
    a_target = fa_target(psin_target)


    data_save=np.column_stack((psin_target, area_target, volume_target, a_target, ne_target, p_target, Te_target, Ti_target, Jtor_target, q_target, psi_target, Ip_target))
    np.savetxt(postproc_folder + 'JOREK2HELENA_profiles.txt', data_save, delimiter=' ', comments='', header = 'psi_n                    area(m^2)                volume(m^3)              minor_radius(m)          density(m^{-3})          presspre(Pa)             Te(eV)                   Ti(eV)                   Jtor(A/m^2)              q                         psi_p                    Ip(A)')  


    plt.figure()
    plt.plot(np.sqrt(psin_q), -q)
    plt.plot(np.sqrt(psin_target), -q_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'$-q$')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_exprs), ne)
    plt.plot(np.sqrt(psin_target), ne_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'density [m$^{-3}$]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_exprs), Ti)
    plt.plot(np.sqrt(psin_target), Ti_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'$T_i$ [eV]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_exprs), Te)
    plt.plot(np.sqrt(psin_target), Te_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'$T_e$ [eV]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_exprs), p)
    plt.plot(np.sqrt(psin_target), p_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'$p$ [Pa]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_exprs), Jtor)
    plt.plot(np.sqrt(psin_target), Jtor_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'$J_\varphi$ [A/m$^2$]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_exprs), psi)
    plt.plot(np.sqrt(psin_target), psi_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'$\psi$ (poloidal magnetic flux)')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_intg3D), volume)
    plt.plot(np.sqrt(psin_target), volume_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'Plasma volume [m$^3$]')
    plt.grid()
    plt.tight_layout()


    plt.figure()
    plt.plot(np.sqrt(psin_intg3D), area)
    plt.plot(np.sqrt(psin_target), area_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'Plasma area [m$^3$]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_intg3D), Ip)
    plt.plot(np.sqrt(psin_target), Ip_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'Plasma current [A]')
    plt.grid()
    plt.tight_layout()

    plt.figure()
    plt.plot(np.sqrt(psin_RZ), a_minor)
    plt.plot(np.sqrt(psin_target), a_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'Minor radius (Rmax - Rmin)/2 [m]')
    plt.grid()
    plt.tight_layout()


    # Read the file
    with open('jorek_36663_3.75-3.95_RZ_Boundary_ed22.txt', 'r') as file:
        lines = file.readlines()

    # Extract numerical values from lines
    R_boundary = np.array([float(val) for val in lines[0].split('=')[1].split(',')])
    Z_boundary = np.array([float(val) for val in lines[1].split('=')[1].split(',')])

    # Read the file
    iq = 8
    formatted_qnum = f"{q_boundary[iq]:.4f}" 
    formatted_psinum = f"{psi_n_boundary[iq]:.5f}" 
    column_names, data_flux = pltm.read_ascii_file(postproc_folder + 'q_surface_' + formatted_qnum + '_s' + str_num + '.txt', 0, droprows = -1)
    R_flux = data_flux[:, 0]
    Z_flux = data_flux[:, 1]
    # theta = math.atan2(R_flux - R_axis, Z_flux - Z_axis)
    angles_rad = np.arctan2(Z_flux - Z_axis, R_flux - R_axis)
    sorted_indices = np.argsort(angles_rad)
    R_flux_sorted = R_flux[sorted_indices]
    Z_flux_sorted = Z_flux[sorted_indices]
    angles_rad_uniform = np.linspace(min(angles_rad), max(angles_rad), len(R_boundary))
    Rflux_target = interp1d(angles_rad[sorted_indices], R_flux_sorted, kind='linear', fill_value="extrapolate")
    R_flux_target = Rflux_target(angles_rad_uniform)
    Zflux_target = interp1d(angles_rad[sorted_indices], Z_flux_sorted, kind='linear', fill_value="extrapolate")
    Z_flux_target = Zflux_target(angles_rad_uniform)


    plt.figure()
    plt.plot(R_boundary, Z_boundary, marker = 'o', linestyle = 'none')
    # plt.plot(R_flux_sorted, Z_flux_sorted, marker = '+', linestyle = 'none')
    plt.plot(R_flux_target, Z_flux_target, marker = '*', linestyle = 'none')
    plt.xlabel(r'R [m]')
    plt.ylabel(r'Z [m]')
    plt.grid()
    plt.tight_layout()

    # Combine into two-column format
    data_boundary = np.column_stack((R_boundary, Z_boundary))
    output_file = postproc_folder + 'JOREK2HELENA_boundary.txt'
    np.savetxt(output_file, data_boundary, header="R_boundary(m)            Z_boundary(m)", delimiter=" ", comments='', )

    data_flux = np.column_stack((R_flux_target, Z_flux_target))
    output_file = postproc_folder + 'JOREK2HELENA_flux_psin_eq_' + formatted_psinum + ' .txt'
    np.savetxt(output_file, data_flux, header="R_flux(m)                Z_flux(m)", delimiter=" ", comments='', )

    data0D = np.column_stack((Iptot, psi_axis, psi_boundary))
    output_file = postproc_folder + 'JOREK2HELENA_0Ddata.txt'
    np.savetxt(output_file, data0D, header="Iptotal(A)               Psi_axis                 Psi_boundary", delimiter=" ", comments='', )

    plt.show()
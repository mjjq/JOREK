#! /usr/bin/env python3

# Read files from jorek2_postproc and generate inputs for HELENA

import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltm
from scipy.interpolate import interp1d

postproc_folder = 'postproc_si_unit/'
postproc_qfile = 'qprofile_s000300.dat'
postproc_exprsfile = 'exprs_averaged_s000300.dat'
postproc_intg3Dfile = 'int3Dprofile_s000300.dat'
postproc_RZfile = 'RZ_maxmin_profile_s000300.dat'
e_charge = 1.602176634e-19 #C
Ip =  8.129874792E+05 # A
psi_boundary = 0.00004
psi_axis = 0.20921

if __name__ == "__main__":

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

    column_names_intg3D, data_intg3D = pltm.read_ascii_file(postproc_folder + postproc_intg3Dfile, 1, droprows = -1)
    psin_intg3D = data_intg3D[:, 0]
    volume = data_intg3D[:, 1]
    area = data_intg3D[:, 2]

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

    fvolume_target = interp1d(psin_intg3D, volume, kind='linear', fill_value="extrapolate")
    volume_target = fvolume_target(psin_target)

    farea_target = interp1d(psin_intg3D, area, kind='linear', fill_value="extrapolate")
    area_target = farea_target(psin_target)

    fa_target = interp1d(psin_RZ, a_minor, kind='linear', fill_value="extrapolate")
    a_target = fa_target(psin_target)


    data_save=np.column_stack((psin_target, area_target, volume_target, a_target, ne_target, p_target, Te_target, Ti_target, Jtor_target, q_target))
    np.savetxt(postproc_folder + 'JOREK2HELENA_profiles.txt', data_save, delimiter=' ', comments='', header = 'psi_n                    area(m^2)                volume(m^3)              minor_radius(m)          density(m^{-3})          presspre(Pa)             Te(eV)                   Ti(eV)                   Jtor(A/m^2)              q')  


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
    plt.plot(np.sqrt(psin_RZ), a_minor)
    plt.plot(np.sqrt(psin_target), a_target, linestyle=':')
    plt.xlabel(r'$\sqrt{\psi_n}$')
    plt.ylabel(r'Minor radius (Rmax - Rmin)/2 [m]')
    plt.grid()
    plt.tight_layout()


    # Read the file
    with open(postproc_folder + 'jorek_36663_3.75-3.95_RZ_Boundary_ed22.txt', 'r') as file:
        lines = file.readlines()

    # Extract numerical values from lines
    R_boundary = np.array([float(val) for val in lines[0].split('=')[1].split(',')])
    Z_boundary = np.array([float(val) for val in lines[1].split('=')[1].split(',')])

    plt.figure()
    plt.plot(R_boundary, Z_boundary, marker = 'o', linestyle = 'none')
    plt.xlabel(r'R [m]')
    plt.ylabel(r'Z [m]')
    plt.grid()
    plt.tight_layout()

    # Combine into two-column format
    data_boundary = np.column_stack((R_boundary, Z_boundary))
    output_file = postproc_folder + 'JOREK2HELENA_boundary.txt'
    np.savetxt(output_file, data_boundary, header="R_boundary(m)            Z_boundary(m)", delimiter=" ", comments='', )

    data0D = np.column_stack((Ip, psi_axis, psi_boundary))
    output_file = postproc_folder + 'JOREK2HELENA_0Ddata.txt'
    np.savetxt(output_file, data0D, header="Ip(A)                    Psi_axis                 Psi_boundary", delimiter=" ", comments='', )

    plt.show()
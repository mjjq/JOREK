#! /usr/bin/env python3
# calculate the ion-diamagnetic frequency at each radial location
# please output the averaged profiles from jorek2_postproc at SI-unit

import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as pltm

m_proton = 1.67262192595e-27 # kg
e_unit   = 1.60217663e-19 # coulomb
r_mode   = 0.3
u0 = 4*np.pi*1e-7         # H/m
ni = 0.98e20              # m^3
rho_mass_axis = ni * m_proton  # kg/m^3
time2si = np.sqrt(u0 * rho_mass_axis) # tau_A => s
omega2si = 1/time2si
print(f'time2si = {time2si}')
B0 = 2.57                 #T

column_names, data = pltm.read_ascii_file('postproc/exprs_averaged_s000000.dat', 1, 0)
rho_p = data[:,0]        # minor radius
ne    = data[:,1]        # m^3
Ti    = data[:,2]        # keV
Te    = data[:,3]        # keV
ptot  = data[:,4]        # Pa
pi    = ne * e_unit * Ti # Pa
dpi_dr1D = np.gradient(pi, rho_p)
idx = np.abs(rho_p - r_mode).argmin()
dpi_dr = dpi_dr1D[idx]
# dpi_dr = 0.5*(1.432-1.183)*1e5/0.2 # Pa/m
rho_mass1D = ne * m_proton
rho_mass = rho_mass1D[idx]

w_star1D = m_proton / e_unit / (rho_mass1D * B0 * rho_p) * dpi_dr1D

w_star = m_proton / e_unit / (rho_mass * B0 * r_mode) * dpi_dr
print(f'w_star = {w_star} [rad/s]')
w_star_norm = w_star / omega2si
print(f'w_star_norm = {w_star_norm} [$\omega_A$]')

print(f'dpi_dr = {dpi_dr} [Pa/m]')
# plt.plot(rho_p, pi)
# plt.plot(rho_p, dpi_dr1D)
plt.plot(rho_p, -w_star1D)
plt.yscale('log')

plt.grid()
plt.xlabel(r'$\sqrt{\psi_n}$')
plt.show()
#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as plm


file_name = 'T_e_modes_n001'

normT = 1.0
for n in range(1, 2):
    for m in range(3):

        start_line = 1 + m * 603
        end_line   = start_line + 601

        column_names, data_mn = plm.read_ascii_file(file_name, start_line, end_line)

        psi_n = data_mn[:,0]
        T     = data_mn[:,1]
        rho_n = np.sqrt(psi_n)

        if(n == 1 and m == 1):
            normT = max(T)

        # plt.plot(rho_n, np.gradient(T, rho_n))
        plt.plot(rho_n, T/normT, label='m/n = '+ str(m) + '/' + str(n))

plt.legend()
plt.xlabel(r'$\sqrt{\psi_n}$')
plt.ylabel(r'$\delta Te_{m,n}$ (A.U.)')
plt.grid()
plt.xlim(0,0.95)
plt.show()
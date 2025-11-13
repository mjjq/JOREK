#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

tauIC = 1.0
rhop = np.linspace(0, 1, 1000)
distance_bnd = 1.0 - rhop
tauIC_shift = 0.25
tauIC_steep = 0.10

tau_IC = tauIC * (0.5e0 - 0.5e0 * np.tanh(-(distance_bnd - tauIC_shift)/tauIC_steep) )

plt.plot(rhop, tau_IC)
plt.xlabel(r'$\rho_p$')
plt.ylabel(r'$\tau_{IC}$ (A.U.)')
plt.show()
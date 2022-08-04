"""
Script for plotting the output from diagnostics/poincare.f90
"""
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt

# Read data
with h5py.File("poincare.h5", "r") as f:
    r      = f["r"][:]
    z      = f["z"][:]
    phi    = f["phi"][:]
    psi    = f["psi"][:]
    iprt   = f["iprt"][:]
    iplane = f["pncrid"][:]
    dist   = f["mileage"][:]

fig = plt.figure()

s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)

# Assign each marker a color when plotting to help in separating field lines
# from one another. These colors are cycled from a group of six colors.
colors = matplotlib.cm.get_cmap('viridis')(np.linspace(0,1,6))

nprt = np.unique(iprt).size
for i in range(nprt):
    idx = np.logical_and.reduce([iprt == i, iplane == 2])
    s1.scatter(r[idx], z[idx], 1, color=colors[np.mod(i,6),:])

for i in range(nprt):
    idx = np.logical_and.reduce([iprt == i, iplane == 1])
    s2.scatter(psi[idx], phi[idx], 1, color=colors[np.mod(i,6),:])

s1.set_xlabel("R [m]")
s1.set_ylabel("z [m]")

s2.set_xlabel(r"$\psi_n$")
s2.set_ylabel("Toroidal angle [rad]")
s2.set_ylim(0,2*np.pi)

plt.show()

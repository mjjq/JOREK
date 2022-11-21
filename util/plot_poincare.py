"""
Script for plotting the output from diagnostics/poincare.f90
"""
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt

# Use color to display the connection length (otherwise it just separates field lines)
color_is_clength = True

# Read data
with h5py.File("poincare.h5", "r") as f:
    r      = f["r"][:]
    z      = f["z"][:]
    phi    = f["phi"][:]
    psi    = f["psi"][:]
    iprt   = f["iprt"][:]
    iplane = f["pncrid"][:]
    dist   = f["mileage"][:]
    mil    = f["mil"][:]

fig = plt.figure()

s1 = fig.add_subplot(1,2,1)
s2 = fig.add_subplot(1,2,2)

if color_is_clength:

    # mileage is the connection length from the marker *initial* position whereas
    # mil stores the distance the marker had travelled up to that point. Therefore
    # we get the connection length at each Poincare crossing as clength = dist - mil
    prt = np.unique(iprt) - 1
    for i in prt:
        idx = iprt == i + 1
        mil[idx] = (dist[i] - mil[idx])

        # If mil reaches mileage, then this marker was confined so set connection length to 0 (appears as NaN after log)
        if np.amin(mil[idx]) == 0: mil[idx] = 0

    cmin = np.nanmin(mil[mil > 0])
    cmax = np.nanmax(mil)
    cmap = mpl.cm.get_cmap('Blues').copy()
    cmap.set_bad("black")

    idx = iplane == 2
    s1.scatter(r[idx], z[idx], 2, mil[idx], norm=mpl.colors.LogNorm(vmin=cmin,vmax=cmax), cmap=cmap)

    idx = iplane == 1
    h2 = s2.scatter(r[idx], phi[idx], 2, mil[idx], norm=mpl.colors.LogNorm(vmin=cmin,vmax=cmax), cmap=cmap)

    cax = plt.colorbar(h2, ax=s2, location='right', extend='both')
    cax.set_label("Connection length [m]")
    
else:
    # Assign each marker a color when plotting to help in separating field lines
    # from one another. These colors are cycled from a group of six colors.
    colors = mpl.cm.get_cmap('viridis')(np.linspace(0,1,6))

    prt = np.unique(iprt) - 1
    for i in prt:
        idx = np.logical_and.reduce([iprt == i, iplane == 2])
        s1.scatter(r[idx], z[idx], 1, color=colors[np.mod(i,6),:])

    for i in prt:
        idx = np.logical_and.reduce([iprt == i, iplane == 1])
        s2.scatter(psi[idx], phi[idx], 1, color=colors[np.mod(i,6),:])

s1.set_xlabel("R [m]")
s1.set_ylabel("z [m]")

s2.set_xlabel(r"$\psi_n$")
s2.set_ylabel("Toroidal angle [rad]")
s2.set_ylim(0,2*np.pi)

plt.show()

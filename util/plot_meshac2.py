#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def fgaus(zs, bgf, xr1, xr2, sig1, sig2, fact):
    znorm1 = 0.39894 / sig1
    znorm2 = 0.39894 / sig2
    zex1 = -0.5 * (zs - xr1) ** 2 / sig1 ** 2
    zex2 = -0.5 * (zs - xr2) ** 2 / sig2 ** 2
    dex1 = -(zs - xr1) / sig1 ** 2
    dex2 = -(zs - xr2) / sig2 ** 2

    f1 = znorm1 * np.exp(zex1)
    f2 = znorm2 * np.exp(zex2)
    df1 = znorm1 * dex1 * np.exp(zex1)
    df2 = znorm2 * dex2 * np.exp(zex2)

    fgaus_val = bgf + (1.0 - bgf) * (f1 + fact * f2) / fact
    dfgauss = (1.0 - bgf) * (df1 + fact * df2) / fact

    return fgaus_val, dfgauss

def meshac2(nr, sg, xr1, xr2, sig1, sig2, bgf, fact):
    if nr <= 0:
        return

    nmax = 10001

    if sg[0] < 0.0 or sg[0] >= 1.0:
        sg[0] = 0.0

    dsi = (1.0 - sg[0]) / (nmax - 1)
    s1 = np.zeros(nmax)
    fsum = np.zeros(nmax)
    s1[0] = sg[0]
    fsum[0] = 0.0
    summ = 0.0
    fint2, dfg = fgaus(s1[0], bgf, xr1, xr2, sig1, sig2, fact)
    
    for i in range(1, nmax):
        s1[i] = sg[0] + i * dsi
        fint1 = fint2
        fint2, dfg = fgaus(s1[i], bgf, xr1, xr2, sig1, sig2, fact)
        summ += (fint1 + fint2) / 2.0 * dsi
        fsum[i] = summ

    fsum /= fsum[-1]
    fsum[-1] = 1.0

    j = 1
    for i in range(1, nr - 1):
        fi = i / (nr - 1)
        while fi > fsum[j] and j < nmax - 1:
            j += 1
        sg[i] = s1[j - 1] + (fi - fsum[j - 1]) / (fsum[j] - fsum[j - 1]) * (s1[j] - s1[j - 1])

    sg[nr - 1] = 1.0
    sg[0] = 0.0


if __name__ == '__main__':
    nradial = 100
    points = np.linspace(1, nradial, nradial)
    radial = np.zeros(nradial)

    meshac2(nradial, radial, xr1 = 0.0, xr2 = 0.2, sig1 = 0.01, sig2 = 0.01, bgf = 0.6, fact = 1.0)
    plt.plot(points, radial, '.', linestyle = 'none')
    plt.xlabel(r'points')
    plt.ylabel(r'$\sqrt{\psi}$')
    plt.grid()
    plt.show()
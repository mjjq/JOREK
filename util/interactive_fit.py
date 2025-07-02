#! /usr/bin/env python3

#
# Purpose: Adjust 1D profile interactively from ASCII files.
#
# Date: 2025-07-02
# Author: Haowei Zhang, IPP Garching
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import sys

# --- Step 0: Parse arguments ---
if len(sys.argv) < 3:
    print("Usage:")
    print("  python interactive_fit.py data.txt --sine")
    print("  python interactive_fit.py data.txt --poly N")
    print("  python interactive_fit.py data.txt --smooth winlen polyorder")
    sys.exit(1)

filename = sys.argv[1]
mode = sys.argv[2]
x_power = float(sys.argv[-1])

# parse options
poly_degree = None
winlen = None
polyorder = None

if mode == '--poly':
    if len(sys.argv) < 4:
        print("Error: You must specify the degree for --poly")
        sys.exit(1)
    poly_degree = int(sys.argv[3])

elif mode == '--smooth':
    if len(sys.argv) < 5:
        print("Error: You must specify window length and polyorder for --smooth")
        sys.exit(1)
    winlen = int(sys.argv[3])
    polyorder = int(sys.argv[4])
elif mode != '--sine':
    print("Error: Unknown mode. Use --sine, --poly N, or --smooth winlen polyorder")
    sys.exit(1)

# --- Step 1: Load data ---
data = np.loadtxt(filename)
x = data[:, 0]**x_power
y = data[:, 1]

delta_y = max(y) - min(y)
delta_x = max(x) - min(x)
accurucy = 0.01 * delta_x
influ_rad = 0.10 * delta_x

# --- Step 2: Define fitting functions ---
def sine_model(x, a, b, c):
    return a * np.sin(b * x + c)

def poly_model(x, *coeffs):
    return np.polyval(coeffs, x)

# --- Step 3: Interactive draggable points ---
class DraggablePoints:
    def __init__(self, ax, x, y, influence_radius=0.5):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.x = x
        self.y = y
        
        self.points, = ax.plot(x, y, 'ro')
        self.points.set_pickradius(5)
        self.points_original, = ax.plot(x, y, 'k--')

        self.selected_marker, = ax.plot([], [], 'yo', markersize=12, zorder=10)

        self.dragging_point = None
        self.initial_y = None
        self.influence_radius = influence_radius

        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        for i, (xi, yi) in enumerate(zip(self.x, self.y)):
            if max(abs(xi - event.xdata), abs(yi - event.ydata)) < accurucy:
                self.dragging_point = i
                self.initial_y = self.y.copy()  # save original y positions
                self.selected_marker.set_data([xi], [yi])
                self.canvas.draw()
                break

    def on_motion(self, event):
        if self.dragging_point is None or event.inaxes != self.ax:
            return

        dy = event.ydata - self.y[self.dragging_point]  # how much y has changed
        x0 = self.x[self.dragging_point]

        # Gaussian influence based on distance from selected point
        influence = np.exp(-((self.x - x0)**2) / (2 * self.influence_radius**2))

        # Apply displacement scaled by influence
        # self.y[self.dragging_point] = self.initial_y[self.dragging_point] + dy * influence[self.dragging_point]
        self.y = self.initial_y + dy * influence

        # Update plots
        self.points.set_ydata(self.y)
        self.selected_marker.set_data([x0], [self.y[self.dragging_point]])
        self.canvas.draw()

    def on_release(self, event):
        self.dragging_point = None
        self.selected_marker.set_data([], [])
        self.canvas.draw()

# --- Step 4: Show draggable plot ---
fig, ax = plt.subplots()
ax.set_title("Drag red points, then close window to continue")
ax.set_xlim(np.min(x) - delta_x * 0.05, np.max(x) + delta_x * 0.05)
ax.set_ylim(np.min(y) - delta_y * 0.50, np.max(y) + delta_y * 0.50)
dr = DraggablePoints(ax, x.copy(), y.copy(), influence_radius= influ_rad)

plt.grid(True)
plt.show()

x_mdf = dr.x
y_mdf = dr.y
x_fit = np.linspace(min(x_mdf), max(x_mdf), 300)

# --- Step 5: Perform smoothing or fitting ---
try:
    if mode == '--sine':
        popt, _ = curve_fit(sine_model, x_mdf, y_mdf, p0=[1, 1, 0], maxfev=10000)
        y_fit = sine_model(x_fit, *popt)
        label = f"Sine Fit: a={popt[0]:.2f}, b={popt[1]:.2f}, c={popt[2]:.2f}"
        print(f"Sine parameters:\na = {popt[0]:.4f}, b = {popt[1]:.4f}, c = {popt[2]:.4f}")

    elif mode == '--poly':
        coeffs = np.polyfit(x_mdf, y_mdf, poly_degree)
        y_fit = np.polyval(coeffs, x_fit)
        coeff_str = ', '.join(f'{c:.3g}' for c in coeffs)
        label = f"Poly{poly_degree} Fit: {coeff_str}"
        print(f"Polynomial coefficients (degree {poly_degree}):")
        print(coeffs)

    elif mode == '--smooth':
        if winlen % 2 == 0:
            winlen += 1  # Savitzky-Golay requires odd window length
            print(f"Window length adjusted to odd number: {winlen}")
        if winlen > len(y_mdf):
            winlen = len(y_mdf) - (1 - len(y_mdf) % 2)
            print(f"Window length reduced to {winlen} to match data length.")

        y_fit = savgol_filter(y_mdf, winlen, polyorder)
        x_fit = x_mdf  # no interpolation
        label = f"Savitzky-Golay Smooth (win={winlen}, order={polyorder})"
        print(f"Applied Savitzky-Golay smoothing with window={winlen}, order={polyorder}")

    # --- Step 6: Plot result ---
    plt.figure()
    plt.plot(x_mdf, y_mdf, 'ro', label='Modified Data')
    plt.plot(x_fit, y_fit, 'b-', label=label)
    plt.title("Curve Result After Adjusting Points")
    plt.legend().set_draggable(True)
    plt.xlim(np.min(x_mdf), np.max(x_mdf))
    plt.ylim(np.min(y_mdf), np.max(y_mdf))
    plt.grid(True)
    plt.show()

    file_save = filename[:-4] + '_mdf.txt'
    print(f'file_save = {file_save}')
    data_save=np.column_stack((x_mdf**(1.0/x_power), y_mdf))
    np.savetxt(file_save, data_save, delimiter=' ', comments='')
    print(f'modified profile of {filename} is saved to {file_save}')

    file_save = filename[:-4] + '_fit.txt'
    print(f'file_save = {file_save}')
    data_save=np.column_stack((x_fit**(1.0/x_power), y_fit))
    np.savetxt(file_save, data_save, delimiter=' ', comments='')
    print(f'modified profile of {filename} is saved to {file_save}')

except RuntimeError as e:
    print(f"Error: {e}")

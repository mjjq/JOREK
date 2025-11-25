#! /usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# -----------------------------
# Read psiRZ_eqdsk
# -----------------------------
def read_psiRZ(filename):
    xx_list, yy_list, psi_list = [], [], []
    with open(filename, 'r') as f:
        lines = f.readlines()

    block = []
    for line in lines:
        line = line.strip()
        if line == "":
            if block:
                xx_list.append(block[0][0])
                yy_list.append([row[1] for row in block])
                psi_list.append([row[2] for row in block])
                block = []
            continue
        block.append(list(map(float, line.split())))

    xx = np.array(xx_list)
    yy = np.array(yy_list)[0]
    psi = np.array(psi_list)
    return xx, yy, psi

# -----------------------------
# Read boundary.dat
# -----------------------------
def read_boundary(filename):
    rb, zb, psi_b = [], [], []
    with open(filename, 'r') as f:
        for line in f:
            vals = list(map(float, line.split()))
            rb.append(vals[0])
            zb.append(vals[1])
            psi_b.append(vals[2])
    return np.array(rb), np.array(zb), np.array(psi_b)

# -----------------------------
# Transpose boundary.dat
# -----------------------------
def save_transposed_boundary(infile="boundary.dat", outfile="R_Z_psi_boundary.dat"):
    # Read the boundary file
    rb, zb, psi_b = [], [], []
    with open(infile, 'r') as f:
        for line in f:
            vals = list(map(float, line.split()))
            rb.append(vals[0])
            zb.append(vals[1])
            psi_b.append(vals[2])

    rb = np.array(rb)
    zb = np.array(zb)
    psi_b = np.array(psi_b)

    # Write transposed format
    with open(outfile, 'w') as f:
        # Convert to comma-separated with same precision
        rb_str = ", ".join(f"{x:.8f}" for x in rb)
        zb_str = ", ".join(f"{x:.8f}" for x in zb)
        psi_str = ", ".join(f"{x:.8f}" for x in psi_b)

        f.write("R_boundary      =  " + rb_str + "\n")
        f.write("Z_boundary      =  " + zb_str + "\n")
        f.write("Psi_boundary    =  " + psi_str + "\n")

    print(f"Saved transposed boundary to {outfile}")
    print(f"n_boundary = {len(rb)}")

# -----------------------------
# Split Te and Ti
# -----------------------------
def split_temperature(infile="jorek_temperature",
                      outfile_Te="jorek_temperature_Te",
                      outfile_Ti="jorek_temperature_Ti",
                      factor=0.5):
    """
    Read psi_ext and T_ext from infile,
    split into Te and Ti using a factor (<1),
    and save to two new files with the same formatting.
    """

    # --------------------------
    # Read the input file
    # --------------------------
    psi_list = []
    T_list = []

    with open(infile, 'r') as f:
        for line in f:
            if line.strip() == "":
                continue
            vals = list(map(float, line.split()))
            psi_list.append(vals[0])
            T_list.append(vals[1])

    psi_ext = np.array(psi_list)
    T_ext = np.array(T_list)

    # --------------------------
    # Split into Te and Ti
    # --------------------------
    Te_ext = T_ext * factor
    Ti_ext = T_ext * (1 - factor)

    # --------------------------
    # Save Te file
    # --------------------------
    with open(outfile_Te, "w") as f:
        for psi, Te in zip(psi_ext, Te_ext):
            f.write(f"{psi:.8e}  {Te:.8e}\n")

    # --------------------------
    # Save Ti file
    # --------------------------
    with open(outfile_Ti, "w") as f:
        for psi, Ti in zip(psi_ext, Ti_ext):
            f.write(f"{psi:.8e}  {Ti:.8e}\n")

    print(f"Saved electron temperature to: {outfile_Te}")
    print(f"Saved ion temperature to:      {outfile_Ti}")

# -----------------------------
# Save q profile
# -----------------------------
def extract_q_profile(infile="profiles_eqdsk", outfile="q_profile.dat_target"):
    """
    Read profiles_eqdsk, take the first and last column,
    and write them into q_profile.dat_target.
    The first line of the file is the data name.
    """
    data = []

    # --- read file ---
    with open(infile, "r") as f:
        for line in f:
            if line.strip() == "":
                continue
            vals = list(map(float, line.split()))
            data.append(vals)

    data = np.array(data)
    psi = data[:, 0]          # first column
    q   = data[:, -1]         # last column

    # --- write new file ---
    with open(outfile, "w") as f:
        f.write("Psi_N   q\n")     # header

        # write psi and q values
        for p, qq in zip(psi, q):
            f.write(f"{p:.8e}  {qq:.8e}\n")

    print(f"Saved q-profile to: {outfile}")

# -----------------------------
# Main
# -----------------------------
xx, yy, psi = read_psiRZ("psiRZ_eqdsk")
rbnd, zbnd, psi_b = read_boundary("boundary.dat")

# Compute ψ_boundary mean and error
psi_b_mean = psi_b.mean()
psi_b_std  = psi_b.std()
print(f"ψ_boundary_mean = {psi_b_mean:.6f}")
print(f"ψ_boundary_std  = {psi_b_std:.2e}")

# Original mesh coordinates
R_orig, Z_orig = np.meshgrid(xx, yy, indexing='ij')

# Flatten for interpolation
points = np.column_stack((R_orig.ravel(), Z_orig.ravel()))
values = psi.ravel()

# Create finer mesh
nR, nZ = 4*len(xx), 4*len(yy)  # 4x resolution increase
R_fine = np.linspace(xx.min(), xx.max(), nR)
Z_fine = np.linspace(yy.min(), yy.max(), nZ)
R_new, Z_new = np.meshgrid(R_fine, Z_fine, indexing='ij')

# Interpolate ψ onto finer mesh
psi_new = griddata(points, values, (R_new, Z_new), method='cubic')

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(8, 10))

# Filled contour
cf = plt.contourf(R_new, Z_new, psi_new, 80, cmap='viridis')
plt.colorbar(cf, label="ψ")

# Contour lines with labels
contours = plt.contour(R_new, Z_new, psi_new, levels=10, colors='k', linewidths=0.8)
plt.clabel(contours, inline=True, fontsize=8, fmt="ψ=%.3f")

# Plot plasma boundary
plt.plot(rbnd, zbnd, 'white', linewidth=2.5)
plt.plot(rbnd, zbnd, 'black', linewidth=1.0)

# Annotate ψ_boundary
txt = f"ψ_boundary = {psi_b_mean:.6f} ± {psi_b_std:.2e}"
plt.text(0.02, 0.02, txt, transform=plt.gca().transAxes,
         fontsize=10, color='yellow', weight='bold',
         bbox=dict(facecolor='black', alpha=0.6, pad=5))

plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.title("High-Resolution ψ(R,Z) with Boundary")
plt.axis("equal")
plt.tight_layout()
plt.show()

save_transposed_boundary("boundary.dat", "R_Z_psi_boundary.dat")

split_temperature(infile="jorek_temperature", 
                  outfile_Te="jorek_temperature_Te",
                  outfile_Ti="jorek_temperature_Ti",
                  factor=0.5)

extract_q_profile(infile="profiles_eqdsk", outfile="q_profile.dat_target")
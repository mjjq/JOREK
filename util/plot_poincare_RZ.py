import numpy as np
import matplotlib.pyplot as plt

# Read data from the file
file_path = 'poinc_R-Z.dat'
data = np.loadtxt(file_path, skiprows=1)  # Skip the first row

# Extract columns R and Z
R = data[:, 0]
Z = data[:, 1]

# Create scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(R, Z, s=4, c='red', marker='.', alpha=0.7)
plt.xlabel('R')
plt.ylabel('Z')
plt.title('Poincaré Plot of R vs Z')
plt.grid(True)

# plt.axis('equal')  # Set the axis scales equally

# Set x and y axis limits based on the data range
plt.xlim(min(R), max(R))
plt.ylim(min(Z), max(Z))

# Adjust the aspect ratio to make axis scales equal
plt.gca().set_aspect('equal', adjustable='box')

plt.show()

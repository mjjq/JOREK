#! /usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the first line to get column names
with open("q_profile.dat", "r") as file:
    first_line = file.readline().strip()

column_names = first_line.split("|")

# Read the data from the file
data = pd.read_csv("q_profile.dat", sep="\s+", skiprows=1, header=None, names=column_names)
x_data = np.sqrt(data[column_names[0]])  # Square root of the first column
y_data = data[column_names[1]]           # Second column

# Plot the data
plt.plot(x_data, y_data, marker='o')
plt.xlabel("Square Root of " + column_names[0])  # X-axis label
plt.ylabel(column_names[1])                      # Y-axis label
plt.title(column_names[1])
plt.grid(True)  # Adding gridlines
# Set x-axis limits
plt.xlim(0, 1)
plt.show()
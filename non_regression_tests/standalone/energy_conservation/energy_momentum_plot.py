#!/usr/bin/env python
from __future__ import print_function

#import h5py
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import glob
import string
import re
import sys


def read_traces(basename):
  # read files of form basename[0-9]*.dat into list of numpy arrays
  files = sorted([f for f in glob.glob('%s*[.]dat'%basename) if re.match(r'^%s\d+[.]dat$'%basename, f)])
  # Filter and sort these files
  return ([np.fromfile(file,sep=' ') for file in files], [int(re.search(r'\d+', s).group()) for s in files])
         

def normalize_traces(traces):
  # Divide all values by the values in the first array to obtain relative values
  # Remove all lost particles
  traces = ma.masked_values(traces,0.0,atol=1e-35)
  return np.divide(traces,traces[0]) # element-wise division


if __name__ == "__main__":
  # Read the type from cli arguments
  basename = sys.argv[1]

  # Read all traces
  traces, t = read_traces(basename)

  # Normalize with the values at t=0
  traces = normalize_traces(traces)

  # Plot overlapping traces with transparency for a measure of the distribution
  plt.plot(t, traces, lw=4, alpha=0.05) # TODO adjust alpha to number of lines
  plt.xlabel("$t$")
  plt.ylabel("increase in %s"%basename);

  # Plot the mean, min and max sharply

  # Save the figure
  plt.savefig('%s_py.png'%basename)

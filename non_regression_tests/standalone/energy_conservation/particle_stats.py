#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import glob
import string
import re
import sys


"""
read files of form basename[0-9]*_[0-9]*.dat into list of numpy arrays
"""
def read_traces(basename):
        filenames = sorted([f for f in glob.glob('%s*_*[.]dat'%basename) if re.match(r'^%s\d+_\d+[.]dat$'%basename, f)])
        # Find number of timesteps by counting number of files produced by cpu 0
        timesteps = [int(re.search(r'\d+', f).group()) for f in filenames if re.match(r'^%s\d+_[0]+[.]dat$'%basename, f)]
        num_steps = len(timesteps)
        # Find number of cpus
        num_cpus = int(len(filenames)/num_steps)

        # Files are sorted, so timesteps are major factor and per-core files are found together
        #  -------> t (axis 1)
        # |X O O O
        # |X O O O
        # |X O O O
        # |O O O O
        # |O O O O
        # |O O O O
        # N (axis 0)
        # Load blocks of one timestep, n_particles/n_cpu and stack them above eachother.
        # Then, stack these side to side to create our matrix.
        return (np.stack([np.concatenate([np.loadtxt(filenames[j]) for j in range(i,i+num_cpus)]) for i in range(0,num_steps)], axis=-1), timesteps)


if __name__ == "__main__":
        # Read the type from cli arguments
        basename = sys.argv[1]

        # Read all traces
        traces, t = read_traces(basename)

        # Mask all zeroes, these are lost particles
        traces = ma.masked_values(traces,0.0,atol=1e-35)

        # Write output files for meta-computations
        np.savetxt("%s_stats.dat.gz"%basename, np.c_[ma.average(traces,axis=1),ma.std(traces,axis=1),traces.min(axis=1),traces.max(axis=1),traces[:,-1],traces[:,0]],
                   header="mean, std, min, max, end, begin")
        print("Wrote %s_stats.dat.gz"%basename)

        # Plot histogram of min/max ratio
        n, bins, patches = plt.hist(traces.min(axis=1)/traces.max(axis=1), bins=20, alpha=0.75)
        plt.xlabel('%s min/max'%basename)
        plt.ylabel('fraction')
        plt.grid(True)
        plt.savefig('%s_minmax.png'%basename)
        print("Wrote %s_minmax.png"%basename)

        # Plot histogram of growth rate
        n, bins, patches = plt.hist((traces[:,-1]-traces[:,0])/ma.average(traces,axis=1), bins=20, alpha=0.75)
        plt.xlabel('%s growth rate'%basename)
        plt.ylabel('fraction')
        plt.grid(True)
        plt.savefig('%s_growth.png'%basename)
        print("Wrote %s_growth.png"%basename)

        # Plot histogram of variance
        n, bins, patches = plt.hist(ma.std(traces,axis=1), bins=20, alpha=0.75)
        plt.xlabel('%s variance'%basename)
        plt.ylabel('fraction')
        plt.grid(True)
        plt.savefig('%s_variance.png'%basename)
        print("Wrote %s_variance.png"%basename)

        # Plot overlapping traces with transparency for a measure of the distribution
        plt.plot(traces.T, lw=4, alpha=0.05)
        plt.xlabel("$t$")
        plt.ylabel("%s"%basename);
        plt.savefig('%s_trace.png'%basename)
        print("Wrote %s_trace.png"%basename)

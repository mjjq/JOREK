#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use("Agg") # Headless backend
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
        return (np.stack([np.concatenate([np.fromfile(filenames[j],dtype='float64') for j in range(i,i+num_cpus)]) for i in range(0,num_steps)], axis=-1), timesteps)


if __name__ == "__main__":
        # Read the type from cli arguments
        basename = sys.argv[1]

        # Read all t
        t, time = read_traces(basename)

        # Mask all zeroes, these are lost particles
        t = ma.masked_values(t,0.0,atol=1e-35)

        # Write output files for meta-computations
        np.savetxt("%s_stats.txt.gz"%basename, np.c_[ma.average(t,axis=1),ma.std(t,axis=1),t.min(axis=1),t.max(axis=1),t[:,0],t[:,-1]],
                   header="mean, std, min, max, begin, end")
        print("Wrote %s_stats.txt.gz"%basename)

        # Write output files with gathered results
        # This does not work with mixed-sign results! (min/max)
        gr = np.fabs(t[:,-1]/t[:,0]-1)
        np.savetxt("%s_results.txt"%basename, np.c_[ma.average(gr),ma.std(gr),gr.min(),gr.max(),np.fabs(ma.average(ma.std(t,axis=1)/ma.average(t,axis=1))),(1-np.fabs(t).min(axis=1)/np.fabs(t).max(axis=1)).max()],
                   header="<end/begin-1>,sd(end/begin-1),min(end/begin-1),max(end/begin-1),sd(t)/mean(t),max(1-min(abs)/max(abs))")
        print("Wrote %s_results.txt"%basename)

        # Plot histogram of min/max ratio
        n, bins, patches = plt.hist(t.min(axis=1)/t.max(axis=1), bins=20, alpha=0.75)
        plt.xlabel('%s min/max'%basename)
        plt.ylabel('fraction')
        plt.grid(True)
        plt.savefig('%s_minmax.png'%basename)
        print("Wrote %s_minmax.png"%basename)
        plt.clf()

        # Plot histogram of growth rate
        n, bins, patches = plt.hist((t[:,-1]-t[:,0])/ma.average(t,axis=1), bins=20, alpha=0.75)
        plt.xlabel('%s growth rate'%basename)
        plt.ylabel('fraction')
        plt.grid(True)
        plt.savefig('%s_growth.png'%basename)
        print("Wrote %s_growth.png"%basename)
        plt.clf()

        # Plot histogram of variance
        n, bins, patches = plt.hist(ma.std(t,axis=1), bins=20, alpha=0.75)
        plt.xlabel('%s stddev'%basename)
        plt.ylabel('fraction')
        plt.grid(True)
        plt.savefig('%s_std.png'%basename)
        print("Wrote %s_std.png"%basename)
        plt.clf()

        # Plot overlapping t with transparency for a measure of the distribution
        fig = plt.plot(time, t.T/t[:,0].T, lw=1, alpha=0.1)
        plt.xlabel("$t$")
        plt.ylabel("%s"%basename);
        plt.savefig('%s_trace.png'%basename)
        print("Wrote %s_trace.png"%basename)
        plt.clf()

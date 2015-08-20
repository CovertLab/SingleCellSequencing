"""
Analysis!

Compute histograms/make scatter plots of the ratio of spike in reads to the number of mapped reads for each cell
"""

"""
Import python packages
"""

import HTSeq 
import collections
import itertools
import os
import subprocess
import collections
import datetime
import yaml
import fnmatch
import shlex
import numpy
import scipy
import scipy.io as sio 
import pyensembl
import h5py
import pandas as pd
import numpy as np
import matplotlib
import cPickle as pickle
matplotlib.use("Agg")
import matplotlib.pyplot as plt

"""
Load all the cells
"""

matplotlib.style.use('ggplot')


direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

spikeins = ['Spike1', 'Spike4', 'Spike7']

spike1_list = []
spike4_list = []
spike7_list = []
total_counts_list = []

for cell in all_cells:
	spike1_list += [cell.spikeins.loc['Spike1']['est_counts']]
	spike4_list += [cell.spikeins.loc['Spike4']['est_counts']]
	spike7_list += [cell.spikeins.loc['Spike7']['est_counts']]	
	total_counts_list += [cell.total_estimated_counts]

spike1_list = np.array(spike1_list)/2
spike4_list = np.array(spike4_list)/2
spike7_list = np.array(spike7_list)/2
total_counts_list = np.array(total_counts_list)/2

sum_spike_list = spike1_list + spike4_list + spike7_list
fraction_spike = sum_spike_list / total_counts_list

fig = plt.figure()
ax = fig.add_subplot(111)

ax.hist(fraction_spike, bins = 40)
ax.set_xlabel('Fraction of estimated counts belonging to RNA spike ins')
ax.set_ylabel('Number of cells')
ax.set_title('Fraction of estimated counts belonging to RNA spike ins')

fig.tight_layout()

plt.savefig("qc3_spikein_fraction_histogram.pdf")

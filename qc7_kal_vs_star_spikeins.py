"""
Analysis!

Compute histograms/make scatter plots of the number of reads for each spike in for each cell
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
from scipy import stats
matplotlib.use("Agg")
import matplotlib.pyplot as plt

"""
Load all the cells
"""

matplotlib.style.use('ggplot')

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_rsem.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

spike1_list = []
spike4_list = []
spike7_list = []

spike1_list_rsem = []
spike4_list_rsem = []
spike7_list_rsem = []


spikeins = ['Spike1', 'Spike4', 'Spike7']

for cell in all_cells:

	spike1_list += [cell.spikeins.loc['Spike1']['est_counts']+1]
	spike4_list += [cell.spikeins.loc['Spike4']['est_counts']+1]
	spike7_list += [cell.spikeins.loc['Spike7']['est_counts']+1]

	spike1_list_rsem += [cell.spikeins_rsem.loc['Spike_1']['est_counts']+1]
	spike4_list_rsem += [cell.spikeins_rsem.loc['Spike_4']['est_counts']+1]
	spike7_list_rsem += [cell.spikeins_rsem.loc['Spike_7']['est_counts']+1]

spike_list = spike1_list + spike4_list + spike7_list
spike_list_rsem = spike1_list_rsem + spike4_list_rsem + spike7_list_rsem

spike1_list = np.array(spike1_list)
spike1_list_rsem = np.array(spike1_list_rsem)

spike4_list = np.array(spike4_list)
spike4_list_rsem = np.array(spike4_list_rsem)

spike7_list = np.array(spike7_list)
spike7_list_rsem = np.array(spike7_list_rsem)

slope, intercept, r_value, p_value, std_err = stats.linregress(spike_list, spike_list_rsem)
x = np.linspace(0,200000,1000)
y = slope *x + intercept

fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111)
ax.scatter(spike1_list, spike1_list_rsem, color = 'b', s = .5, alpha = 1, label = 'Spike 1')
ax.scatter(spike4_list, spike4_list_rsem, color = 'k', s = .5, alpha = 1, label = 'Spike 4')
ax.scatter(spike7_list, spike7_list_rsem, color = 'g', s = .5, alpha = 1, label = 'Spike 7')


ax.plot(x,y,'r')
ax.text(.2,.9, 'r = ' + str(r_value), ha='center', va='center', transform=ax.transAxes)
ax.set_xlim([.8,1e6])
ax.set_ylim([.8,1e6])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Kallisto counts + 1', fontsize = 12)
ax.set_ylabel('STAR/RSEM counts + 1', fontsize = 12)
ax.set_title('Kallisto vs STAR/RSEM spike in counts', y= 1.05, fontsize = 14)

plt.legend(loc = 4)

fig.tight_layout()
plt.savefig("plots/qc7_kalvsstar_spikeins.pdf")



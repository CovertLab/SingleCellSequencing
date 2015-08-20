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
matplotlib.use("Agg")
import matplotlib.pyplot as plt

"""
Load all the cells
"""

matplotlib.style.use('ggplot')

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

spike1_list = []
spike4_list = []
spike7_list = []

spikeins = ['Spike1', 'Spike4', 'Spike7']

for cell in all_cells:

	spike1_list += [cell.spikeins.loc['Spike1']['est_counts']]
	spike4_list += [cell.spikeins.loc['Spike4']['est_counts']]
	spike7_list += [cell.spikeins.loc['Spike7']['est_counts']]

spike1_list = np.array(spike1_list)/2
spike4_list = np.array(spike4_list)/2
spike7_list = np.array(spike7_list)/2

fig = plt.figure()
ax = fig.add_subplot(111)

ax.patch.set_facecolor('white')
spk1 = fig.add_subplot(311)
spk4 = fig.add_subplot(312)
spk7 = fig.add_subplot(313)

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor = 'w', top = 'off', bottom = 'off', left = 'off', right = 'off')

spk1.hist(spike1_list, bins = 40, label = 'Spike 1 reads')
spk1.legend()

spk4.hist(spike4_list, bins = 40, label = 'Spike 4 reads', color = 'red')
spk4.legend()

spk7.hist(spike7_list, bins = 40, label = 'Spike 7 reads', color = 'green')
spk7.legend()


ax.set_xlabel('Estimated counts')
ax.set_ylabel('Number of cells')
ax.set_title('Estimated counts for RNA spike ins')


plt.legend()

fig.tight_layout()
plt.savefig("qc2_spikein_histogram.pdf")

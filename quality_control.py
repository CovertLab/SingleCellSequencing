"""
Quality control

Remove low quality cells from the list of cells
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
import seq_functions

matplotlib.style.use('ggplot')

"""
Load cells
"""

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

all_cells_qc = seq_functions.quality_control(all_cells)

"""
Plot histograms of mapped vs unmapped cells for cells post qc filtering
"""

unmapped_list = []
mapped_list = []

for cell in all_cells_qc:
	mapped_list += [cell.num_mapped]
	unmapped_list += [cell.num_unmapped]

mapped = np.array(mapped_list)
unmapped = np.array(unmapped_list)

mp = plt.hist(mapped, bins = 40, label = 'Mapped reads')
unmp = plt.hist(unmapped, bins = 40, label = 'Unmapped reads')
plt.xlabel('Number of reads')
plt.ylabel('Number of cells')
plt.title('Number of reads mapping to the transcriptome')
plt.legend()
plt.savefig("plots/qc1_num_mapped_histogram_postqc.pdf")

"""
Plot spike in histograms for cells post qc filtering
"""

spike1_list = []
spike4_list = []
spike7_list = []

spikeins = ['Spike1', 'Spike4', 'Spike7']

for cell in all_cells_qc:

	spike1_list += [cell.spikeins.loc['Spike1']['est_counts']]
	spike4_list += [cell.spikeins.loc['Spike4']['est_counts']]
	spike7_list += [cell.spikeins.loc['Spike7']['est_counts']]

spike1_list = np.array(spike1_list)
spike4_list = np.array(spike4_list)
spike7_list = np.array(spike7_list)

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
plt.savefig("plots/qc2_spikein_histogram_postqc.pdf")

"""
Plot spike in fraction for cells post qc filtering
"""

spikeins = ['Spike1', 'Spike4', 'Spike7']

spike1_list = []
spike4_list = []
spike7_list = []
total_counts_list = []

for cell in all_cells_qc:
	spike1_list += [cell.spikeins.loc['Spike1']['est_counts']]
	spike4_list += [cell.spikeins.loc['Spike4']['est_counts']]
	spike7_list += [cell.spikeins.loc['Spike7']['est_counts']]	
	total_counts_list += [cell.total_estimated_counts]

spike1_list = np.array(spike1_list)
spike4_list = np.array(spike4_list)
spike7_list = np.array(spike7_list)
total_counts_list = np.array(total_counts_list)

sum_spike_list = spike1_list + spike4_list + spike7_list
fraction_spike = sum_spike_list / total_counts_list

fig = plt.figure()
ax = fig.add_subplot(111)

ax.hist(fraction_spike, bins = 40)
ax.set_xlabel('Fraction of estimated counts belonging to RNA spike ins')
ax.set_ylabel('Number of cells')
ax.set_title('Fraction of estimated counts belonging to RNA spike ins')

fig.tight_layout()

plt.savefig("plots/qc3_spikein_fraction_histogram_postqc.pdf")

"""
Save cells post qc filtering in pickle object
"""

file_name_save = os.path.join(direc, 'all_cells_qc.pkl')

pickle.dump(all_cells_qc, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

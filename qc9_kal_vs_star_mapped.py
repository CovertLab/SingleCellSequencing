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

kal_list = []
rsem_list = []

index_rsem = all_cells[0].fpkm_rsem.index
index_kal = all_cells[0].fpkm.index

joint_index = list(set(index_rsem).intersection(set(index_kal)))
joint_index = joint_index[0:500]
# print joint_index

cell = all_cells[0]
# print list(cell.fpkm_rsem.loc[joint_index] +1)
counter = 0
for cell in all_cells:
	print counter
	counter += 1
	rsem_list += list([cell.num_mapped_rsem]) 
	kal_list += list([cell.num_mapped])

rsem_list = np.array(rsem_list)
kal_list = np.array(kal_list)


slope, intercept, r_value, p_value, std_err = stats.linregress(kal_list, rsem_list)
x = np.linspace(0,1400000,1000)
y = x


fig = plt.figure(figsize = (6,6))
ax = fig.add_subplot(111)
# ax.scatter(rsem_list, spike1_list_rsem, color = 'b', s = .5, alpha = 1, label = 'Spike 1')

# ax.set_xlim([.8,1e5])
# ax.set_ylim([.8,1e5])
ax.scatter(kal_list, rsem_list, color = 'b', s = .5, alpha = 1)
ax.plot(x,y,'r')
# ax.text(.2,.9, 'r = ' + str(r_value), ha='center', va='center', transform=ax.transAxes)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlabel('Kallisto reads', fontsize = 12)
ax.set_ylabel('STAR/RSEM reads', fontsize = 12)
ax.set_title('Kallisto vs STAR/RSEM mapped reads comparison', y= 1.05, fontsize = 14)

ax.set_ylim([0,1400000])
ax.set_xlim([0,1400000])
ax.set_xticks([0,400000,800000,1200000])
ax.set_yticks([0,400000,800000,1200000])



# plt.legend(loc = 4)

fig.tight_layout()
plt.savefig("plots/qc9_kalvsstar_mapped.pdf")



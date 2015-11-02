"""
Analysis!

Plot the frequency of dropout events as a function of average gene expression
"""

"""
Import python packages
"""

import HTSeq 
import time
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
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc.pkl'

all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))
num_of_cells = len(all_cells)

cell = all_cells[0]
gene_keys = cell.tpm.index

dropout_frequency = {}
mean_tpm = {}

dropout_frequency_list = []
mean_tpm_list = []

for gene in gene_keys:
	num_dropout = 0
	total_tpm = 0
	for cell in all_cells:
		if cell.tpm.loc[gene] == 0:
			num_dropout += 1
		total_tpm += cell.tpm.loc[gene]

	dropout_frequency[gene] = num_dropout/num_of_cells
	mean_tpm[gene] = total_tpm/num_of_cells
	dropout_frequency_list += [num_dropout/num_of_cells]
	mean_tpm_list += [total_tpm/num_of_cells]

"""
Plot dependence of dropout frequency on expression level
"""

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(mean_tpm_list, dropout_frequency_list)

ax.set_xlabel('Mean TPM across all cells')
ax.set_ylabel('Dropout frequency')
ax.set_title('Dropout frequency vs mean expression')

plt.savefig("qc4_dropout_tpm.pdf")

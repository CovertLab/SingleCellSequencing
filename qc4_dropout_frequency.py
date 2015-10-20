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

matplotlib.style.use('ggplot')


direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc.pkl'

all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))
num_of_cells = len(all_cells)

cell = all_cells[0]
gene_keys = cell.tpm.index


dropout_frequency = pd.DataFrame(np.zeros((len(gene_keys),2), dtype = 'float32'), index = gene_keys, columns = ['tpm','dropout'])

for cell in all_cells:
	dropout_frequency.loc[:,'tpm'] += (1+cell.tpm)/num_of_cells
	zero_genes = cell.tpm[cell.tpm == 0].index.tolist()
	dropout_frequency.ix[zero_genes,'dropout'] += np.float32(num_of_cells) ** -1

mean_tpm_list = np.array(dropout_frequency['tpm'].tolist())
dropout_frequency_list = np.array(dropout_frequency['dropout'].tolist())

# print np.where(mean_tpm_list > 100)[0].shape
# print np.where(dropout_frequency_list < 0.1)[0].shape

print sum(dropout_frequency_list)
print 'The number of genes with a TPM+1 > 100 is ' + str(len(np.where(mean_tpm_list > 100)[0]))
print 'The number of genes with a dropout frequency < 0.1 is ' + str(len(np.where(dropout_frequency_list < 0.1)[0]))
# print dropout_frequency

"""
Plot dependence of dropout frequency on expression level
"""

fig = plt.figure(figsize = (5,4))
ax = fig.add_subplot(111)
ax.scatter(mean_tpm_list, dropout_frequency_list, color = 'b', s = .5, alpha = 1)
ax.set_xscale('log')
ax.set_xlabel('Mean (TPM+1) across all cells', fontsize = 12)
ax.set_ylabel('Dropout frequency', fontsize = 12)
ax.set_title('Dropout frequency vs mean expression', y= 1.05, fontsize = 14)
ax.set_ylim([-0.05,1.05])
ax.set_yticks([0,.2,.4,.6,.8,1])
ax.set_xlim([1,1e5])
fig.tight_layout()
plt.savefig("plots/qc4_dropout_tpm.pdf")

"""
Limit transcriptomes to genes with high enough expression so that the dropout effect is small
"""

# dropout_threshold = 0.1

# genes_with_low_dropout = []

# for gene in gene_keys:
# 	if dropout_frequency[gene] < dropout_threshold:
# 		genes_with_low_dropout += [gene]
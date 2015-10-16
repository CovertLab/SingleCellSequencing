"""
Compute the geometric mean size factors to normalize the transcriptome fpkm's
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
all_cell_file = 'all_cells_qc.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Compute geometric mean normalization for each cells - Remove genes with 0 counts prior to computing size factor
"""

cell = all_cells[0]
gene_keys = list(cell.fpkm.index)

"""
Compute reference sample for each gene
"""

reference_sample = {}
for gene in gene_keys:
	running_product = 1
	number_of_cells = 0
	for cell in all_cells:
		if cell.fpkm.loc[gene] > 0:
			running_product *= cell.fpkm.loc[gene]
			number_of_cells += 1
	print gene
	if number_of_cells > 0:
		reference_sample[gene] = running_product ** (1/number_of_cells)
	else:
		reference_sample[gene] = 0

"""
Compute size factor for each cell
"""

for cell in all_cells:
	list_of_genes = []
	for gene in gene_keys:
		if cell.fpkm.loc[gene] > 0:
			if reference_sample[gene] > 0:
				list_of_genes += [cell.fpkm.loc[gene] / reference_sample[gene]]
	cell.size_factor = np.median(list_of_genes)

for cell in all_cells:
	print cell.size_factor

file_name_save = os.path.join(direc, 'all_cells_norm.pkl')

pickle.dump(all_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)


 


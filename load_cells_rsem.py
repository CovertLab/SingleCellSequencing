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
import cPickle as pickle
import seq_functions

direc = '/scratch/PI/mcovert/dvanva/sequencing/'

# Load matfiles
mat_file_path = os.path.join(direc,'imaging/')
matfiles = seq_functions.load_matfiles(mat_file_path)
dynamics_data = seq_functions.dynamics_class(matfiles)

counted_direc = 'counted'
counted_rsem_direc = 'counted_rsem'
librarys_to_load = ['library1', 'library2', 'library3', 'library4', 'library5', 'library6', 'library7','library8', 'library9', 'library10']

all_cells = []
all_cells_justkal = []
for library in librarys_to_load:
	print 'Loading ' + library
	lib_path = os.path.join(direc,library,counted_direc)
	lib_path_rsem = os.path.join(direc,library, counted_rsem_direc)
	file_list = os.listdir(lib_path)
	for h5_file in file_list:
		if fnmatch.fnmatch(h5_file, r'*.h5'):
			print h5_file
			h5_rsem = h5_file[:-3] + '.isoforms.results.h5'
			print h5_rsem

			cell = seq_functions.cell_object_rsem(h5_kallisto_file = os.path.join(lib_path, h5_file), h5_rsem_file = os.path.join(lib_path_rsem, h5_rsem), dictionary = dynamics_data)
			cell2 = seq_functions.cell_object(h5_file = os.path.join(lib_path, h5_file), dictionary = dynamics_data)
			all_cells += [cell]
			all_cells_justkal += [cell2]
			
print len(all_cells)

all_cells_qc = seq_functions.quality_control(all_cells)
all_cells_justkal_qc = seq_functions.quality_control(all_cells_justkal)

print len(all_cells_qc), len(all_cells_justkal_qc)

all_cells_qc = seq_functions.remove_unidentified_genes(all_cells_qc)
all_cells_qc = seq_functions.remove_jackpotting_genes(all_cells_qc)
all_cells_qc = seq_functions.add_tpm_normalization(all_cells_qc)
all_cells_qc = seq_functions.add_tpm_normalization_rsem(all_cells_qc)

all_cells_qc = seq_functions.remove_low_expression_genes(all_cells_qc)


file_name_save = os.path.join(direc, 'all_cells_rsem_qc.pkl')
pickle.dump(all_cells_qc, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)





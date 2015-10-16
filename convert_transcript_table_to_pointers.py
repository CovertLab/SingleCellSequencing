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

# Load list of cell objects

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

num_cells = len(all_cells)
est_counts = all_cells[0].transcripts.loc[:,'est_counts']
fpkm = all_cells[0].transcripts.loc[:,'fpkm']
tpm = all_cells[0].transcripts.loc[:,'tpm']
all_cells[0].transcript_loc = 0

all_cells[0].transcripts = None

for j in xrange(1,num_cells):
	print "Processing cell " + str(j) + " of " + str(num_cells)
	cell = all_cells[j]
	est_counts_list = [est_counts,cell.transcripts.loc[:,'est_counts']]
	fpkm_list = [fpkm,cell.transcripts.loc[:,'fpkm']]
	tpm_list = [tpm,cell.transcripts.loc[:,'tpm']]

	est_counts = pd.concat(est_counts_list)
	fpkm = pd.concat(fpkm_list)
	tpm = pd.concat(tpm_list) 

	all_cells[j].transcripts = None
	all_cells[j].transcript_loc = j

file_name_save = os.path.join(direc, 'all_cells_pointers.pkl')
pickle.dump(all_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

file_name_save = os.path.join(direc, 'est_counts.pkl')
pickle.dump(est_counts, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

file_name_save = os.path.join(direc, 'fpkm.pkl')
pickle.dump(fpkm, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

file_name_save = os.path.join(direc, 'tpm.pkl')
pickle.dump(tpm, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)


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
librarys_to_load = ['library1', 'library2', 'library3', 'library4', 'library5', 'library6', 'library7','library8', 'library9', 'library10']

all_cells = []

for library in librarys_to_load:
	print 'Loading ' + library
	lib_path = os.path.join(direc,library,counted_direc)
	file_list = os.listdir(lib_path)
	for h5_file in file_list:
		if fnmatch.fnmatch(h5_file, r'*.h5'):
			cell = seq_functions.cell_object(h5_file = os.path.join(lib_path, h5_file), dictionary = dynamics_data)
			all_cells += [cell]

file_name_save = os.path.join(direc, 'all_cells.pkl')

pickle.dump(all_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)


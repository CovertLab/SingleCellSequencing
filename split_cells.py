"""
Cell list splitter

Split list of cells by time point to make the pickle objects faster to load

Only select stimulated cells
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

times = [0, 75, 150, 300]

all_cells_stim = []
for cell in all_cells:
	if cell.condition == 'Stim':
		all_cells_stim += [cell]

for t in times:
	all_cells_time = []

	for cell in all_cells_stim:
		if cell.time_point == t:
			all_cells_time += [cell]

	file_name_save = os.path.join(direc, 'all_cells_' + str(t) + 'min.pkl')

	pickle.dump(all_cells_time, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

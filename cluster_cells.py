"""
Cluster time series

Load the 0 time point and another timepoint and cluster the cells by timeseries using 
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



times = [75, 150, 300]

for t in times:
	print t

	"""
	Load cells
	"""
	direc = '/scratch/PI/mcovert/dvanva/sequencing/'
	all_cell_0min_file = 'all_cells_0min.pkl'
	all_cell_time_file = 'all_cells_' + str(t) + 'min.pkl'

	all_cells_0 = pickle.load(open(os.path.join(direc,all_cell_0min_file)))
	all_cells_t = pickle.load(open(os.path.join(direc,all_cell_time_file)))

	all_cells = all_cells_0 + all_cells_t

	all_cells = seq_functions.cell_cluster(all_cells)

	file_name_save = os.path.join(direc, 'all_cells_' + str(t) + 'min_clustered.pkl')

	pickle.dump(all_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)





"""
Analysis!

Plot a histogram of size factors for all of the cells
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

"""
Load all the cells
"""

matplotlib.style.use('ggplot')

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_norm.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Histogram the size factors
"""

size_factors = []

for cell in all_cells:
	size_factors += [cell.size_factor]

sf = plt.hist(size_factors, bins = 40)
plt.xlabel('Size factor')
plt.ylabel('Number of cells')
plt.title('Size factor')
plt.savefig("plots/trial_6_size_factor.pdf")
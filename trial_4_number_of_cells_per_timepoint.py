"""
Analysis!

Make a histogram of the number of cells in each time point
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

matplotlib.style.use('ggplot')
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Figure out the length of the longest time trace
"""

longest_time = 0
number_of_cells = 0
t = [0, 75, 150, 300]
t_label = ['0 min', '75 min', '150 min', '300 min']
counter_list = []

for j in xrange(4):
	for cell in all_cells:
		if cell.time_point == t[j]:
			number_of_cells += 1
	counter_list += [number_of_cells]
	number_of_cells = 0

"""
Plot histogram
"""
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(4)
width = 0.35
ax.set_xticks(ind + width/2)
ax.bar(ind, counter_list, width)
xtickNames = ax.set_xticklabels(t_label)
plt.setp(xtickNames, rotation = 45, fontsize = 10, ha = 'right')
ax.set_ylabel('Number of cells')
ax.set_title('Number of cells per time point - ' + str(len(all_cells)) + ' cells total')

plt.savefig("trial_4_number_of_cells_per_timepoint.pdf")




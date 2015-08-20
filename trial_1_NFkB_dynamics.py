"""
Analysis!

Plot the dynamics for the first 20 cells in the all_cell object
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
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

cell_counter = 0 

fig = plt.figure()
ax = fig.add_subplot(111)

for cell in all_cells:
	if cell_counter > 20:
		break
	ax.plot(cell.NFkB_dynamics, label = 'Cell' + str(cell_counter))
	ax.set_xlabel('Time (timepoints)')
	ax.set_ylabel('Normalized fluorescence (au)')
	ax.set_title('Nuclear concentration of NFkB')
	cell_counter += 1


plt.savefig("trial_1_NFkBdynamics_20cells.pdf")




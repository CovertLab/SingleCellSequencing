"""
Analysis!

Plot a heatmap for all of the cells in the all cells list
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

# matplotlib.style.use('ggplot')
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

def cleanAxis(ax):
	ax.set_frame_on(False)
	for label in ax.axes.get_xticklabels():
		label.set_visible(False)
	for label in ax.axes.get_yticklabels():
		label.set_visible(False)
	for tick in ax.axes.get_xticklines():
		tick.set_visible(False)
	for tick in ax.axes.get_yticklines():
		tick.set_visible(False)
	for spine in ax.spines.values():
		spine.set_visible(False)

"""
Figure out the length of the longest time trace
"""

longest_time = 0
number_of_cells = 0
t = 75

for cell in all_cells:
	if cell.time_point == t:
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])

heat_map = np.zeros((number_of_cells,longest_time))

"""
Fill up the heat map matrix
"""

cell_counter = 0
for cell in all_cells:
	if cell.time_point == t:
		number_of_cells += 1
		dynam = cell.NFkB_dynamics
		heat_map[cell_counter,0:dynam.shape[0]] = dynam
		cell_counter += 1

fig = plt.figure(figsize = (6,8))
ax = fig.add_subplot(111)
cleanAxis(ax)

cax = ax.imshow(heat_map, cmap = plt.get_cmap('Reds'), interpolation = 'none')
ax.set_xlabel('Time')
ax.set_ylabel('Cells')
ax.set_title('75 minute NFkB activity heatmap - ' + str(number_of_cells) + ' cells', y = 1.05)
fig.colorbar(cax, ticks = [0, 1], orientation = 'vertical')

plt.savefig("trial_3_dynamics_heatmap_75min.pdf")





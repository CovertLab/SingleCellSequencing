"""
Analysis!

Cluster the time traces and then plot a heatmap for the dynamics traces
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
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()

# matplotlib.style.use('ggplot')
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_w_jackpot.pkl'
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


times = [75] #[75, 150, 300]

R = rpy2.robjects.r
DTW = importr('dtw')

for t in times:

	longest_time = 0
	number_of_cells = 0
	for cell in all_cells:
		if cell.time_point == t and cell.condition == 'Stim':
			number_of_cells += 1
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])

	dynamics_matrix = np.zeros((number_of_cells,longest_time))

	"""
	Fill up the heat map matrix
	"""

	cell_counter = 0
	for cell in all_cells:
		if cell.time_point == t and cell.condition == 'Stim':
			dynam = cell.NFkB_dynamics
			dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
			cell_counter += 1

	"""
	Perform hierarchical clustering
	"""

	# distance_matrix = np.zeros((number_of_cells, number_of_cells))
	# for i in xrange(number_of_cells):
	# 	print i
	# 	for j in xrange(number_of_cells):
	# 		# print i, j
	# 		# distance_matrix[i,j] = np.linalg.norm(dynamics_matrix[i,:] - dynamics_matrix[j,:])
	# 		alignment = R.dtw(dynamics_matrix[i,:], dynamics_matrix[j,:], keep = True)
	# 		distance_matrix[i,j] = alignment.rx('distance')[0][0]

	# np.savez('/home/dvanva/SingleCellSequencing/150_dynamics_distance_matrix.npz', distance_matrix = distance_matrix)
	dynamics_load = np.load('/home/dvanva/SingleCellSequencing/75_dynamics_distance_matrix.npz')
	distance_matrix = dynamics_load['dynamics_distance_matrix']
	Y = sch.linkage(distance_matrix, method = 'ward')

	"""
	Plot dendrogram
	"""

	fig = plt.figure()
	ax_dendro = fig.add_axes([0.09, 0.1, 0.2, 0.8], frame_on = False)
	Z = sch.dendrogram(Y, orientation = 'right', color_threshold = 0.3*np.amax(Y[:,2]))

	ax_dendro.set_xticks([])
	ax_dendro.set_yticks([])

	"""
	Plot heatmap
	"""

	ax_heatmap = fig.add_axes([0.3, 0.1, 0.6, 0.8])
	index = Z['leaves']
	dynamics_ordered = dynamics_matrix[index,:]
	im = ax_heatmap.matshow(dynamics_ordered, aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')
	fig.colorbar(im, ticks = [0, 1], orientation = 'vertical')

	ax_heatmap.set_title(str(t) +' minute NFkB activity heatmap - ' + str(number_of_cells) + ' cells', y = 1.05)
	ax_heatmap.set_xlabel('Time')
	ax_heatmap.set_yticks([])
	ax_heatmap.set_xticks([])

	plt.savefig("plots/trial_5_dynamics_clustering_" + str(t) + "min.pdf")





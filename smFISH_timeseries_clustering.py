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

direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name_save = os.path.join(direc, "good_cells_300min.pkl")
all_cells = pickle.load(open(file_name_save))

times = [300]

R = rpy2.robjects.r
DTW = importr('dtw')
DTWCLUST = importr('dtwclust')

for t in times:
	longest_time = 0
	number_of_cells = 0
	for cell in all_cells:
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.norm_med.shape[0]])

	dynamics_matrix = np.zeros((number_of_cells,longest_time))

	"""
	Fill up the heat map matrix
	"""

	cell_counter = 0
	for cell in all_cells:
		dynam = cell.norm_med
		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		cell_counter += 1

	"""
	Perform hierarchical clustering
	"""

	distance_matrix = np.zeros((number_of_cells, number_of_cells))
	for i in xrange(number_of_cells):
		print i
		for j in xrange(number_of_cells):
			alignment = R.SBD(dynamics_matrix[i,:], dynamics_matrix[j,:], znorm = True)
			distance_matrix[i,j] = alignment.rx('dist')[0][0]

	np.savez('/scratch/PI/mcovert/dvanva/sequencing/smFISH/'+str(t)+"_dynamics_distance_matrix_kshape.npz", distance_matrix = distance_matrix)
	dynamics_load = np.load('/scratch/PI/mcovert/dvanva/sequencing/smFISH/'+str(t)+"_dynamics_distance_matrix_kshape.npz")
	distance_matrix = dynamics_load['distance_matrix']
	Y = sch.linkage(distance_matrix, method = 'ward')
	ind_dynamics = sch.fcluster(Y,0.5*np.amax(Y[:,2]),'distance') - 1


	"""
	Plot dendrogram
	"""

	fig = plt.figure()
	ax_dendro = fig.add_axes([0.09, 0.1, 0.2, 0.8], frame_on = False)
	Z = sch.dendrogram(Y, orientation = 'right', color_threshold = 0.5*np.amax(Y[:,2]))

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

	plt.savefig("plots/smFISH_dynamics_clustering_" + str(t) + "min.pdf")

	"""
	Assign clusters
	"""
	for j in xrange(number_of_cells):
		all_cells[j].clusterID = ind_dynamics[j]

	file_name_save = os.path.join(direc, "good_cells_" + str(t) + "min_cluster_w_smFISH_traces.pkl")
	pickle.dump(all_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)


	"""
	Compute k-shape average for each cluster
	"""
	max_clust_id = np.amax(ind_dynamics)
	min_clust_id = np.amin(ind_dynamics)
	cluster_id_list = np.arange(min_clust_id, max_clust_id + 1)

	cluster_dynamics_DBA = np.zeros((max_clust_id-min_clust_id+1,longest_time), dtype = 'float32')

	for j in xrange(len(cluster_id_list)):
		cluster_id = cluster_id_list[j]
		num_of_cells = len(np.where(ind_dynamics == cluster_id)[0])
		dyn_matrix = np.zeros((num_of_cells,longest_time))
		
		dyn_list = []

		cell_counter = 0
		for cell in all_cells:
			if cell.clusterID == cluster_id:
				dyn_matrix[cell_counter,0:dynam.shape[0]] = cell.norm_med
				dyn_list += [cell.norm_med]
				cell_counter += 1
		# temp = dba(dyn_list)
		temp = R.shape_extraction(dyn_matrix, znorm = True)
		temp = np.array(temp)
		temp -= temp[0]

		temp /= np.amax(temp) 
		cluster_dynamics_DBA[j,:] = temp
		# print cluster_dynamics_DBA

	colors = ['g', 'r', 'c','purple','yellow']
	np.savez("/scratch/PI/mcovert/dvanva/sequencing/smFISH/300_cluster_avg_kshape_smFISH.npz", cluster_dynamics_avg = cluster_dynamics_DBA)


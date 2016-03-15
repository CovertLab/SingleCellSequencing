"""
smFISH data analysis
Perform analysis for single cell analysis data
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
# import h5py
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
from seq_functions import smFISH_cell
import rpy2
from rpy2.robjects.packages import importr
import cPickle as pickle
from seq_functions import smFISH_cell_minimal
rpy2.robjects.numpy2ri.activate()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42

"""
Initialize R instances
"""

R = rpy2.robjects.r
DTW = importr('dtw')
DTWCLUST = importr('dtwclust')

"""
Load excel files
"""

direc = "/scratch/PI/mcovert/dvanva/sequencing/cell_density"

list_of_dfs = []
file_name = os.path.join(direc, "LPS2500_df.csv")
list_of_dfs += [pd.read_csv(file_name)]
file_name = os.path.join(direc, "LPS5000_df.csv")
list_of_dfs += [pd.read_csv(file_name)]
file_name = os.path.join(direc, "LPS10000_df.csv")
list_of_dfs += [pd.read_csv(file_name)]
file_name = os.path.join(direc, "LPS15000_df.csv")
list_of_dfs += [pd.read_csv(file_name)]
file_name = os.path.join(direc, "LPS25000_df.csv")
list_of_dfs += [pd.read_csv(file_name)]

list_of_densities = ["2500", "5000", "10000", "15000", "25000"]

"""
Load cluster averages
"""
direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name = os.path.join(direc,"300_cluster_avg_kshape_c1.npz")
file_load = np.load(file_name)
cluster_dynamics_avg = file_load["cluster_dynamics_avg"]



dense_count = 0
for df in list_of_dfs:
	list_of_cells = []
	for row in xrange(len(df.index)):
		norm_med = df.iloc[row,1:].values
		cell = smFISH_cell_minimal(norm_med = norm_med, cluster_dynamics_avg = cluster_dynamics_avg)
		list_of_cells += [cell]

	"""
	Plot heat map using the c1 clustering
	"""

	"""
	Fill up the heat map matrix
	"""

	longest_time = 0
	number_of_cells = 0
	for cell in list_of_cells:
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.norm_med.shape[0]])

	dynamics_matrix = np.zeros((number_of_cells,longest_time))

	cluster_len_1 = 0
	cluster_len_2 = 0
	cluster_len_3 = 0

	for cell in list_of_cells:
		if cell.clusterID == 0:
			cluster_len_1 += 1
		if cell.clusterID == 1:
			cluster_len_2 += 1
		if cell.clusterID == 2:
			cluster_len_3 += 1

	print cluster_len_1, cluster_len_2, cluster_len_3
	frac_1 = np.float(cluster_len_1) / np.float(cluster_len_1 + cluster_len_2 + cluster_len_3)
	frac_2 = np.float(cluster_len_2) / np.float(cluster_len_1 + cluster_len_2 + cluster_len_3)
	frac_3 = np.float(cluster_len_3) / np.float(cluster_len_1 + cluster_len_2 + cluster_len_3)

	print frac_1, frac_2, frac_3

	cell_counter = 0

	fig = plt.figure()
	ax_heatmap_1 = fig.add_axes([0.3, 0.1 +0.005, 0.6, 0.8*frac_1 - 0.005], frame_on = True)
	ax_heatmap_2 = fig.add_axes([0.3, 0.1 +0.005+ 0.8*frac_1, 0.6 , 0.8*frac_2 - 0.005], frame_on = True)
	ax_heatmap_3 = fig.add_axes([0.3, 0.1 +0.005+ 0.8*frac_1 + 0.8*frac_2, 0.6, 0.8*frac_3 - 0.005], frame_on = True)

	for cell in list_of_cells:
		if cell.clusterID == 0:
			dynam = cell.norm_med
			dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
			cell_counter += 1

	for cell in list_of_cells:
		if cell.clusterID == 1:
			dynam = cell.norm_med
			dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
			cell_counter += 1

	for cell in list_of_cells:
		if cell.clusterID == 2:
			dynam = cell.norm_med
			dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
			cell_counter += 1


	# ax_heatmap = fig.add_axes([0.3, 0.1, 0.6, 0.8])
	im1 = ax_heatmap_1.matshow(dynamics_matrix[0:cluster_len_1,:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')
	im2 = ax_heatmap_2.matshow(dynamics_matrix[cluster_len_1:cluster_len_1+cluster_len_2,:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')
	im3 = ax_heatmap_3.matshow(dynamics_matrix[cluster_len_1+cluster_len_2:,:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')


	# fig.colorbar(im1, ticks = [0, 1], orientation = 'vertical')

	# ax_heatmap.xaxis.set_ticks(np.arange(63))
	ax_heatmap_1.xaxis.set_ticklabels([])
	ax_heatmap_1.xaxis.set_ticks([])
	ax_heatmap_1.yaxis.set_ticklabels([])
	ax_heatmap_1.yaxis.set_ticks([])
	ax_heatmap_2.yaxis.set_ticklabels([])
	ax_heatmap_2.yaxis.set_ticks([])
	ax_heatmap_2.xaxis.set_ticklabels([])
	ax_heatmap_2.xaxis.set_ticks([])
	ax_heatmap_3.yaxis.set_ticklabels([])
	ax_heatmap_3.yaxis.set_ticks([])
	ax_heatmap_3.xaxis.set_ticklabels([])
	ax_heatmap_3.xaxis.set_ticks([])
	ax_heatmap_3.set_title("Cell density - " + list_of_densities[dense_count] + " - 300 minutes - " + str(number_of_cells) + ' cells', y = 1.05)
	ax_heatmap_1.set_xlabel('Time (minutes)')

	plt.savefig("plots/300min_cell_density_dynamics_c1_clustering_"+ list_of_densities[dense_count] +".pdf")
	dense_count += 1





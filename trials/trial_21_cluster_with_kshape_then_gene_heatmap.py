"""
Analysis!

Cluster the time traces and then compare the gene expression for each cluster
"""

"""
Import python packages
"""

import HTSeq 
import time
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
import matplotlib as mpl
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2 import robjects as ro

from dba import dba
from dba import align_to

from rpy2.robjects.vectors import DataFrame as RDataFrame
from rpy2 import rinterface
from rpy2.robjects import conversion

@conversion.py2ro.register(pd.DataFrame)
def py2ro_pandasdataframe(obj):
    ri_dataf = conversion.py2ri(obj)
    # cast down to an R list (goes through a different code path
    # in the DataFrame constructor, avoiding `str(k)`) 
    ri_list = rinterface.SexpVector(ri_dataf)
    return RDataFrame(ri_list)

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
# mpl.style.use('ggplot')

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')
scde = importr("scde")

# Load data sets in R
r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")
r("""load("/scratch/PI/mcovert/dvanva/sequencing/counts_data.RData")""")
r("o.fpm = scde.expression.magnitude(o.ifm, counts = counts_data_int)")
# r("print(o.fpm)")
# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

# Determine which genes to look at
genes_to_plot = ["Prdx1", "Ccl5", "Ccl3", "Saa3", "Il1f9"]

"""
Analyze all the time points
"""

times_to_analyze = [300]
for time_point in times_to_analyze:

	print "Analyzing " + str(time_point) + " minute time point"
	all_cells = []
	cell_names = []
	longest_time = 0
	number_of_cells = 0

	for cell in all_cells_total:
		if cell.time_point == time_point and cell.condition == 'Stim':
			number_of_cells += 1
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]
			cell_names += [cell.id]

	dynamics_matrix = np.zeros((number_of_cells,longest_time), dtype = 'float32')

	"""
	Fill up the dynamics heat map matrix
	"""
	cell_counter = 0
	for cell in all_cells:
		dynam = cell.NFkB_dynamics
		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		cell_counter += 1

	"""
	Perform hierarchical clustering of the dynamics
	"""
	distance_matrix_dynamics = np.zeros((number_of_cells, number_of_cells))

	dynamics_load = np.load('/home/dvanva/SingleCellSequencing/' + str(time_point)+'_dynamics_distance_matrix_kshape.npz')
	distance_matrix_dynamics = dynamics_load["distance_matrix"]

	Y_dynamics = sch.linkage(distance_matrix_dynamics, method = 'ward')
	ind_dynamics = sch.fcluster(Y_dynamics,0.5*np.amax(Y_dynamics[:,2]),'distance')
	cluster_list = np.arange(np.amin(ind_dynamics), np.amax(ind_dynamics)+1)

	for j in xrange(number_of_cells):
		all_cells[j].clusterID = ind_dynamics[j]

	cluster_dict = {}
	cluster_name_dict = {}

	for cell in all_cells:
		cluster_dict[cell.id] = str(cell.clusterID)

	for cluster in cluster_list:
		cluster_name_dict[str(cluster)] = []
		for cell in all_cells:
			if cell.clusterID == cluster:
				cluster_name_dict[str(cluster)] += [cell.id]

	genes_matrix = np.zeros((number_of_cells,len(genes_to_plot)), dtype = 'float32')
	counter1 = 0
	for cell in all_cells:
		counter2 = 0
		for gene in genes_to_plot:
			cell_name = "'" + cell.id + "'"
			gene_name = "'" + gene + "'"
			fpm = np.log2(np.float32(pandas2ri.ri2py(r("exp(o.fpm[" + gene_name + "," + cell_name + "])")))+1)
			genes_matrix[counter1,counter2] = fpm

			counter2 += 1
		counter1 += 1


	"""
	Plot dendrogram
	"""

	fig = plt.figure()
	ax_dendro = fig.add_axes([0.09, 0.1, 0.2, 0.8], frame_on = False)
	Z = sch.dendrogram(Y_dynamics, orientation = 'right', color_threshold = 0.5*np.amax(Y_dynamics[:,2]))

	ax_dendro.set_xticks([])
	ax_dendro.set_yticks([])

	"""
	Plot heatmap
	"""

	cluster_len_1 = np.float32(len(cluster_name_dict["1"]))
	cluster_len_2 = np.float32(len(cluster_name_dict["2"]))
	cluster_len_3 = np.float32(len(cluster_name_dict["3"]))

	frac_1 = cluster_len_1 / (cluster_len_1 + cluster_len_2 + cluster_len_3)
	frac_2 = cluster_len_2 / (cluster_len_1 + cluster_len_2 + cluster_len_3)
	frac_3 = cluster_len_3 / (cluster_len_1 + cluster_len_2 + cluster_len_3)

	print frac_1, frac_2, frac_3
	ax_heatmap_1 = fig.add_axes([0.3, 0.1+0.005, 0.6, 0.8*frac_1 - 0.01], frame_on = True)
	ax_heatmap_2 = fig.add_axes([0.3, 0.1 +0.005+ 0.8*frac_1, 0.6 , 0.8*frac_2 - 0.01], frame_on = True)
	ax_heatmap_3 = fig.add_axes([0.3, 0.1 +0.005+ 0.8*frac_1 + 0.8*frac_2, 0.6, 0.8*frac_3 - 0.01], frame_on = True)

	index = Z['leaves']
	genes_ordered = genes_matrix[index,:]
	im = ax_heatmap_1.matshow(genes_ordered[0:cluster_len_1], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('coolwarm'), interpolation = 'none')
	im = ax_heatmap_2.matshow(genes_ordered[cluster_len_1:cluster_len_2], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('coolwarm'), interpolation = 'none')
	im = ax_heatmap_3.matshow(genes_ordered[cluster_len_2:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('coolwarm'), interpolation = 'none')

	# fig.colorbar(im, ticks = [0, 20], orientation = 'vertical')

	for i in xrange(genes_matrix.shape[1]):
		print genes_to_plot[i]
		ax_heatmap_1.text(i-.1, -6, '' + genes_to_plot[i], rotation = 270)

	ax_heatmap_3.set_title("300 min gene expression heatmap", y = 1.05)
	ax_heatmap_1.set_yticks([])
	ax_heatmap_1.set_xticks([])
	ax_heatmap_2.set_yticks([])
	ax_heatmap_2.set_xticks([])
	ax_heatmap_3.set_yticks([])
	ax_heatmap_3.set_xticks([])

	plt.savefig("plots/trial_21_300min.pdf")


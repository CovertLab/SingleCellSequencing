"""
Analysis!

Cluster the time traces and then plot the DBA for each cluster traces
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
import matplotlib
import cPickle as pickle
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
from dba import dba, k_dba, align_to
from dtw import dtw as _dtw
from functools import partial

_RELATIVE_TOLERANCE = 1e-2
_MAX_ITERS = 15
_WARP_PENALTY = 2

dtw = partial(_dtw, warp_penalty = _WARP_PENALTY)

mpl.rcParams['pdf.fonttype'] = 42
# matplotlib.style.use('ggplot')
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_w_jackpot.pkl'

all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))


"""
Figure out the length of the longest time trace
"""

time_point = 300
n_clust = 4


all_cells = []
longest_time = 0
number_of_cells = 0
for cell in all_cells_total:
	if cell.time_point == time_point and cell.condition == 'Stim':
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
		all_cells += [cell]
dynamics_matrix = np.zeros((number_of_cells,longest_time), dtype = 'float32')


"""
Fill up the dynamics heat map matrix
"""
start_time = time.time()

cell_counter = 0
dyn_list = []
for cell in all_cells:
	dynam = cell.NFkB_dynamics
	dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
	dyn_list += [dynam]
	cell_counter += 1

time_diff = str(round(time.time()-start_time,1))
print 'Dynamics matrix filled in %s seconds' % time_diff


"""
Perform DBA/kmeans clustering of the dynamics
"""
start_time = time.time()
centers, errors, weights = k_dba(dyn_list, n_clust)

time_diff = str(round(time.time()-start_time,1))
print 'Dynamics clustering completed in %s seconds' % time_diff

# print centers
# print len(errors[0]), len(weights[0]), len(all_cells)
# print weights

# Assign cluster ID
for j in xrange(number_of_cells):
	clusters = np.zeros((n_clust,))
	for i in xrange(n_clust):
		clusters[i] = weights[i][j]
	all_cells[j].clusterID = np.argmax(clusters)

"""
Plot each cluster and associated DBA average
"""

colors = ['g', 'r', 'b', 'k']
cluster_id_list = np.arange(0,n_clust)
for j in xrange(len(cluster_id_list)):
	
	times = np.arange(0,longest_time)*5 #[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]

	plt.clf()

	for cell in all_cells:
		if cell.clusterID == cluster_id_list[j]:
			# plt.plot(times,cell.NFkB_dynamics, color = str(.8))
			plt.plot(times,align_to(centers[j],cell.NFkB_dynamics), color = str(.8))

	plt.plot(times, centers[j], color = colors[j], linewidth = 2, label = 'Cluster ' + str(j+1))

	plt.xlabel('Time (min)', fontsize = 16)
	plt.ylabel('Nuclear localization (au)', fontsize = 16)
	plt.xlim([0, times[-1]])
	plt.xticks([0,times[-1]],  fontsize = 16)
	plt.ylim([0, 1.05])
	plt.yticks([0, 1],  fontsize = 16)

	plt.tight_layout()
	file_name = "DBA_kmeans_average_cluster_" + str(time_point) + "min_" + str(j) + ".pdf"
	plt.savefig("plots/DBA/" + file_name)

plt.clf()
for j in xrange(len(cluster_id_list)):
	
	times = np.arange(0,longest_time)*5 #[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]


	plt.plot(times, centers[j], color = colors[j], linewidth = 2, label = 'Cluster ' + str(j+1))

	plt.xlabel('Time (min)', fontsize = 16)
	plt.ylabel('Nuclear localization (au)', fontsize = 16)
	plt.xlim([0, times[-1]])
	plt.xticks([0,times[-1]],  fontsize = 16)
	plt.ylim([0, 1.05])
	plt.yticks([0, 1],  fontsize = 16)

plt.tight_layout()
file_name = "DBA_kmeans_average_cluster_" + str(time_point) + "min" + ".pdf"
plt.savefig("plots/DBA/" + file_name)

	






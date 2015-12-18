"""
Analysis!

Cluster the time traces and then plot the k-shape average for each cluster traces
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
from dba import dba
from dba import align_to

mpl.rcParams['pdf.fonttype'] = 42
# matplotlib.style.use('ggplot')
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'

all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))


"""
Figure out the length of the longest time trace
"""

time_point = [300]

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')

all_cells = []
longest_time = 0
number_of_cells = 0
for t in time_point:
	for cell in all_cells_total:
		if cell.time_point == t and cell.condition == 'Stim':
			number_of_cells += 1
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]
dynamics_matrix = np.zeros((number_of_cells,longest_time), dtype = 'float32')

"""
Fill up the dynamics heat map matrix
"""
start_time = time.time()

cell_counter = 0
for cell in all_cells:
	dynam = cell.NFkB_dynamics
	dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
	cell_counter += 1

time_diff = str(round(time.time()-start_time,1))
print 'Dynamics matrix filled in %s seconds' % time_diff

"""
Perform hierarchical clustering of the dynamics
"""
start_time = time.time()
distance_matrix_dynamics = np.zeros((number_of_cells, number_of_cells))

dynamics_load = np.load('/home/dvanva/SingleCellSequencing/' + str(time_point[0])+'_dynamics_distance_matrix_kshape.npz')
if time_point[0] == 75:
	distance_matrix_dynamics = dynamics_load["distance_matrix"]
if time_point[0] == 150:
	distance_matrix_dynamics = dynamics_load["distance_matrix"]
if time_point[0] == 300:
	distance_matrix_dynamics = dynamics_load["distance_matrix"]

Y_dynamics = sch.linkage(distance_matrix_dynamics, method = 'ward')
ind_dynamics = sch.fcluster(Y_dynamics,0.5*np.amax(Y_dynamics[:,2]),'distance')

time_diff = str(round(time.time()-start_time,1))
print 'Dynamics clustering completed in %s seconds' % time_diff

for j in xrange(number_of_cells):
	all_cells[j].clusterID = ind_dynamics[j]

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
			dyn_matrix[cell_counter,0:dynam.shape[0]] = cell.NFkB_dynamics
			dyn_list += [cell.NFkB_dynamics]
			cell_counter += 1
	# temp = dba(dyn_list)
	temp = R.shape_extraction(dyn_matrix, znorm = True)
	temp = np.array(temp)
	temp -= temp[0]

	temp /= np.amax(temp) 
	cluster_dynamics_DBA[j,:] = temp
	# print cluster_dynamics_DBA

colors = ['g', 'r', 'c','purple','yellow']
for j in xrange(len(cluster_id_list)):
	
	times = np.arange(0,longest_time)*5
	plt.clf()

	# for cell in all_cells:
	# 	if cell.clusterID == cluster_id_list[j]:
	# 		plt.plot(times,align_to(cluster_dynamics_DBA[j,:],cell.NFkB_dynamics), color = str(.8))

	plt.plot(times, cluster_dynamics_DBA[j,:], color = colors[j], linewidth = 2, label = 'Cluster ' + str(j+1))

	plt.xlabel('Time', fontsize = 16)
	plt.ylabel('Nuclear localization (au)', fontsize = 16)
	plt.xlim([0, times[-1]])
	plt.xticks([],  fontsize = 16)
	plt.ylim([-0.05,1.05])
	plt.yticks([0,1],  fontsize = 16)
	# plt.ylim([np.amin(cluster_dynamics_DBA)*1.05, np.amax(cluster_dynamics_DBA)*1.05])
	# plt.yticks([np.amin(cluster_dynamics_DBA)*1.05, np.amax(cluster_dynamics_DBA)*1.05],  fontsize = 16)

	plt.tight_layout()
	file_name = "kshape_all_clusters_"+str(time_point[0])+"min_cluster_" + str(j) + ".pdf"
	plt.savefig("plots/kshape/" + file_name)
	






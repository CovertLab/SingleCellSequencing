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

mpl.rcParams['pdf.fonttype'] = 42
# matplotlib.style.use('ggplot')

direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name_save = os.path.join(direc, "good_cells_300min.pkl")
all_cells = pickle.load(open(file_name_save))

"""
Figure out the length of the longest time trace
"""

time_point = [300]

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')

longest_time = 0
number_of_cells = 0
for t in time_point:
	for cell in all_cells:
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.norm_med.shape[0]])

dynamics_matrix = np.zeros((number_of_cells,longest_time), dtype = 'float32')

"""
Fill up the dynamics heat map matrix
"""
start_time = time.time()

cell_counter = 0
for cell in all_cells:
	dynam = cell.norm_med
	dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
	cell_counter += 1

time_diff = str(round(time.time()-start_time,1))
print 'Dynamics matrix filled in %s seconds' % time_diff

"""
Perform hierarchical clustering of the dynamics
"""
start_time = time.time()

dynamics_load = np.load('/scratch/PI/mcovert/dvanva/sequencing/smFISH/'+str(t)+"_dynamics_distance_matrix_kshape.npz")
distance_matrix = dynamics_load['distance_matrix']

Y_dynamics = sch.linkage(distance_matrix, method = 'ward')
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

cluster_dynamics = np.zeros((max_clust_id-min_clust_id+1,longest_time), dtype = 'float32')

for j in xrange(len(cluster_id_list)):
	cluster_id = cluster_id_list[j]
	num_of_cells = len(np.where(ind_dynamics == cluster_id)[0])
	dyn_matrix = np.zeros((num_of_cells,longest_time))
	
	dyn_list = []

	cell_counter = 0
	for cell in all_cells:
		if cell.clusterID == cluster_id:
			dyn_matrix[cell_counter,0:dynam.shape[0]] = cell.norm_med
			cell_counter += 1

	temp = R.shape_extraction(dyn_matrix, znorm = True)
	temp = np.array(temp)
	temp -= temp[0]

	temp /= np.amax(temp) 

	cluster_dynamics[j,:] = temp

colors = ['g', 'r', 'c','purple','yellow']
# np.savez('/scratch/PI/mcovert/dvanva/sequencing/smFISH/'+str(time_point[0])+"_cluster_avg_kshape.npz", cluster_dynamics_avg = cluster_dynamics)

for j in xrange(len(cluster_id_list)):
	
	times = np.arange(0,longest_time)*5
	# plt.clf()

	# for cell in all_cells:
	# 	if cell.clusterID == cluster_id_list[j]:
	# 		plt.plot(times,align_to(cluster_dynamics_DBA[j,:],cell.NFkB_dynamics), color = str(.8))

	plt.plot(times, cluster_dynamics[j,:], color = colors[j], linewidth = 2, label = 'Cluster ' + str(j+1))

	plt.xlabel('Time', fontsize = 16)
	plt.ylabel('Nuclear localization (au)', fontsize = 16)
	plt.xlim([0, times[-1]])
	plt.xticks([],  fontsize = 16)
	plt.ylim([-0.5,1.05])
	plt.yticks([0,1],  fontsize = 16)
	# plt.ylim([np.amin(cluster_dynamics_DBA)*1.05, np.amax(cluster_dynamics_DBA)*1.05])
	# plt.yticks([np.amin(cluster_dynamics_DBA)*1.05, np.amax(cluster_dynamics_DBA)*1.05],  fontsize = 16)

	plt.tight_layout()
file_name = "smFISH_kshape_all_clusters_"+str(time_point[0])+"min_cluster_" + str(j) + ".pdf"

plt.savefig("plots/" + file_name)
	






"""
Analysis!

Cluster the time traces and then plot a heatmap for the dynamics traces
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

# matplotlib.style.use('ggplot')
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'

all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

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


times = [75]

R = rpy2.robjects.r
DTW = importr('dtw')

all_cells = []
longest_time = 0
number_of_cells = 0
for t in times:
	for cell in all_cells_total:
		if cell.time_point == t:
			number_of_cells += 1
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]

zero_times = [0]
all_cells_zero = []
for t in zero_times:
	for cell in all_cells_total:
		if cell.time_point == t:
			all_cells_zero += [cell]
num_of_cells_zero = len(all_cells_zero)

print len(all_cells)

dynamics_matrix = np.zeros((number_of_cells,longest_time))

list_of_genes = all_cells[0].tpm.index.get_values()
number_of_genes = len(list_of_genes)
genes_matrix = np.zeros((number_of_cells, number_of_genes))


zero_tpm_mean = pd.DataFrame(np.zeros((len(list_of_genes)), dtype = 'float32'), index = list_of_genes, columns = ['tpm'])

for cell in all_cells_zero:
	zero_tpm_mean.loc[:,'tpm'] += cell.tpm/num_of_cells_zero


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
for i in xrange(number_of_cells):
	print str(i) + ' of ' + str(number_of_cells)
	for j in xrange(number_of_cells):
		alignment = R.dtw(dynamics_matrix[i,:], dynamics_matrix[j,:], keep = True)
		distance_matrix_dynamics[i,j] = alignment.rx('distance')[0][0]

ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
Y_dynamics = sch.linkage(distance_matrix_dynamics, method = 'centroid')
Z_dynamics = sch.dendrogram(Y_dynamics, orientation = 'right')
ind_dynamics = sch.fcluster(Y_dynamics,3,'maxclust')

ax1.set_xticks([])
ax1.set_yticks([])
time_diff = str(round(time.time()-start_time,1))
print 'Dynamics clustering completed in %s seconds' % time_diff

for j in xrange(number_of_cells):
	all_cells[j].clusterID = ind_dynamics[j]

"""
Regroup cells into clusters
"""
all_clusters = []
for cluster in ind_dynamics:
	temp_cluster = []
	for cell in all_cells:
		if cell.clusterID = cluster:
			temp_cluster += [cell]

	tpm_mean = pd.DataFrame(np.zeros((len(list_of_genes)), dtype = 'float32'), index = list_of_genes, columns = ['tpm'])
	for cell in temp_cluster:
		tpm_mean.loc[:,'tpm'] += cell.tpm/len(temp_cluster)


"""
Fill up the gene heat map matrix
"""
start_time = time.time()

print np.squeeze(np.array(zero_tpm_mean.loc[list_of_genes]+1)).shape
print np.array(cell.tpm.loc[list_of_genes]+1).shape

cell_counter = 0
for cell in all_cells:
	genes_matrix[cell_counter,:] = np.array(cell.tpm.loc[list_of_genes]+1) / np.squeeze(np.array(zero_tpm_mean.loc[list_of_genes]+1))
	cell_counter += 1

genes_matrix = np.log2(genes_matrix)

time_diff = str(round(time.time()-start_time,1))
print 'Gene matrix filled in %s seconds' % time_diff

"""
Set up figure 
"""

cmap=plt.cm.coolwarm

### Scale the max and min colors so that 0 is white/black
vmin=genes_matrix.min()
vmax=genes_matrix.max()
# vmax = max([vmax,abs(vmin)])
# vmin = vmax*-1
norm = mpl.colors.Normalize(vmin, vmax) ### adjust the max and min to scale these colors

### Scale the Matplotlib window size
default_window_hight = 8.5
default_window_width = 12
fig = plt.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
color_bar_w = 0.015 ### Sufficient size to show
	
## calculate positions for all elements
# ax1, placement of dendrogram 1, on the left of the heatmap
[ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
width_between_ax1_axr = 0.004
height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix

# axr, placement of row side colorbar
[axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
axr_x = ax1_x + ax1_w + width_between_ax1_axr
axr_y = ax1_y; axr_h = ax1_h
width_between_axr_axm = 0.004

# axc, placement of column side colorbar
[axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
axc_x = axr_x + axr_w + width_between_axr_axm
axc_y = ax1_y + ax1_h + height_between_ax1_axc
height_between_axc_ax2 = 0.004

# axm, placement of heatmap for the data matrix
[axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
axm_x = axr_x + axr_w + width_between_axr_axm
axm_y = ax1_y; axm_h = ax1_h
axm_w = axc_w

# ax2, placement of dendrogram 2, on the top of the heatmap
[ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
ax2_x = axr_x + axr_w + width_between_axr_axm
ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
ax2_w = axc_w

# axcb - placement of the color legend
[axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.09]


"""
Perform hierarchical clustering of the transcriptome
"""
start_time = time.time()
distance_matrix_genes = np.zeros((number_of_genes, number_of_genes))

for i in xrange(len(list_of_genes)):
	print str(i) + ' of ' + str(len(list_of_genes))
	for j in xrange(len(list_of_genes)):
		a = genes_matrix[:,i]
		b = genes_matrix[:,j]
		distance_matrix_genes[i,j] = np.linalg.norm(a-b)

ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
Y_genes = sch.linkage(distance_matrix_genes, method = 'centroid')
Z_genes = sch.dendrogram(Y_genes, color_threshold = 0.4*max(Y_genes[:,2]))
ind_genes = sch.fcluster(Y_genes, 0.2*max(Y_genes[:,2]),'distance')
ax2.set_xticks([])
ax2.set_yticks([])
time_diff = str(round(time.time()-start_time,1))
print 'Gene clustering completed in %s seconds' % time_diff


# Plot distance matrix.
axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
data_matrix = genes_matrix

idx2 = Z_genes['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
data_matrix = data_matrix[:,idx2]
ind_genes = ind_genes[idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram

idx1 = Z_dynamics['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
data_matrix = data_matrix[idx1,:]   # xt is transformed x
ind_dynamics = ind_dynamics[idx1] ### reorder the flat cluster to match the order of the leaves the dendrogram

print ind_dynamics

### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
im = axm.matshow(data_matrix, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
axm.set_xticks([]) ### Hides x-ticks
axm.set_yticks([])

# Plot colside colors
# axc --> axes for column side colorbar
axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
dc = numpy.array(ind_genes, dtype=int)
dc.shape = (1,len(ind_genes)) 
im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
axc.set_xticks([]) ### Hides ticks
axc.set_yticks([])

# Plot rowside colors
# axr --> axes for row side colorbar
dynamics_ordered = dynamics_matrix[ind_dynamics,:]

axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
dr = numpy.array(ind_dynamics, dtype=int)
dr.shape = (len(ind_dynamics),1)
im_r = axr.matshow(dynamics_ordered, aspect='auto', origin='lower', cmap=plt.get_cmap('Reds'))
axr.set_xticks([]) ### Hides ticks
axr.set_yticks([])

# Plot color legend
axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
axcb.tick_params(labelsize = 8)
axcb.set_title("colorkey")

### Render the graphic
# if len(row_header)>50 or len(column_header)>50:
# 	plt.rcParams['font.size'] = 5
# else:
# 	plt.rcParams['font.size'] = 8

filename = 'plots/trial_7_dual_clustering_75min_foldchange.pdf'
plt.savefig(filename)
print 'Exporting:',filename
plt.savefig(filename) #,dpi=200
plt.show()





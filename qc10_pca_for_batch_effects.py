"""
Analysis!

Perform PCA to look for batch effects
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
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
from sklearn.decomposition import PCA

rpy2.robjects.numpy2ri.activate()
matplotlib.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
# matplotlib.style.use('ggplot')

"""
Load cells
"""
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_w_jackpot.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Create data matrix - each row corresponds to a different cell, each column to the TPM of a different gene
"""

num_of_cells = len(all_cells)
list_of_genes = all_cells[0].tpm.index.get_values()
num_of_genes = len(list_of_genes)

list_of_cell_ids = []
for cell in all_cells:
	list_of_cell_ids += [cell.id]

data_matrix = pd.DataFrame(np.zeros((num_of_cells, num_of_genes), dtype = 'float32'), index = list_of_cell_ids, columns = list_of_genes)
data_matrix_postpc = pd.DataFrame(np.zeros((num_of_cells, 4), dtype = 'float32'), index = list_of_cell_ids, columns = [0,1,2,3])

for cell in all_cells:
	data_matrix.loc[cell.id,:] = np.log2(cell.tpm.transpose()+1)

pca = PCA(n_components = 4)
data_matrix_pc = pca.fit_transform(data_matrix)
data_matrix_postpc.loc[:,:] = data_matrix_pc

print data_matrix_postpc
print pca.explained_variance_ratio_

"""
Plot pca components
"""

comp_0 = data_matrix_postpc.loc[:,0]
comp_1 = data_matrix_postpc.loc[:,1]

fig = plt.figure(figsize = (6,6))
ax = fig.add_subplot(111)
ax.scatter(comp_0, comp_1, color = 'b', s = .5, alpha = 1)

ax.set_xlabel('Component 0', fontsize = 12)
ax.set_ylabel('Component 1', fontsize = 12)
ax.set_title('PCA of Kallisto counts', y= 1.05, fontsize = 14)

plt.savefig("plots/qc10_pca.pdf")

"""
Plot 3D pca components and label by time point
"""

times = [0,75,150,300]
colors_dict = {0:'b',75:'r',150:'g',300:'k'}

blue_line = mpl.lines.Line2D([],[],color = 'blue', label = '0 minutes')
red_line = mpl.lines.Line2D([],[],color = 'red', label = '75 minutes')
green_line = mpl.lines.Line2D([],[],color = 'green', label = '150 minutes')
black_line = mpl.lines.Line2D([],[],color = 'black', label = '300 minutes')

comp_0 = []
comp_1 = []
comp_2 = []
colors = []

for cell in all_cells:
	comp_0 += [data_matrix_postpc.loc[cell.id,0]]
	comp_1 += [data_matrix_postpc.loc[cell.id,1]]
	comp_2 += [data_matrix_postpc.loc[cell.id,2]]
	colors += [colors_dict[cell.time_point]]

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(comp_0, comp_1, comp_2, color = colors, s = 2, alpha = 1)
ax.set_xlim([-80,40])
ax.set_ylim([-30,50])
ax.set_zlim([-30,30])
ax.set_xlabel('Component 0', fontsize = 12)
ax.set_ylabel('Component 1', fontsize = 12)
ax.set_zlabel('Component 2', fontsize = 12)

ax.set_title('PCA of Kallisto counts', y= 1.1, fontsize = 14)
ax.legend(handles = [blue_line, red_line, green_line, black_line], loc = 2)
plt.savefig("plots/qc10_pca_color_by_time.pdf")

"""
Plot pca components and label by library number
"""

times = [0,75,150,300]
colors_dict = {0:'b',75:'r',150:'g',300:'k'}

comp_0 = []
comp_1 = []
label_list = []

for cell in all_cells:
	comp_0 += [data_matrix_postpc.loc[cell.id,0]]
	comp_1 += [data_matrix_postpc.loc[cell.id,1]]
	label_list += [cell.library_number]
comp_0 = data_matrix_postpc.loc[:,0]
comp_1 = data_matrix_postpc.loc[:,1]

fig = plt.figure(figsize = (10,10))
ax = fig.add_subplot(111)
ax.scatter(comp_0, comp_1, color = 'b', s = 0.5, alpha = 1)
for i,j,k in zip(comp_0,comp_1,label_list):
	ax.annotate(str(k), xy = (i,j), fontsize = 4)

ax.set_xlabel('Component 0', fontsize = 12)
ax.set_ylabel('Component 1', fontsize = 12)
ax.set_title('PCA of Kallisto counts', y= 1.05, fontsize = 14)
plt.savefig("plots/qc10_pca_color_by_library_number.pdf")

"""
Plot 3D pca components and label by library number (libraries 1-5)
"""

times = [0,75,150,300]
colors_dict = {1:'b',2:'r',3:'g',4:'k',5:'c'}

blue_line = mpl.lines.Line2D([],[],color = 'blue', label = 'Library 1')
red_line = mpl.lines.Line2D([],[],color = 'red', label = 'Library 2')
green_line = mpl.lines.Line2D([],[],color = 'green', label = 'Library 3')
black_line = mpl.lines.Line2D([],[],color = 'black', label = 'Library 4')
cyan_line = mpl.lines.Line2D([],[],color = 'cyan', label = 'Library 5')


comp_0 = []
comp_1 = []
comp_2 = []
colors = []

for cell in all_cells:
	if cell.library_number in [1,2,3,4,5]:
		comp_0 += [data_matrix_postpc.loc[cell.id,0]]
		comp_1 += [data_matrix_postpc.loc[cell.id,1]]
		comp_2 += [data_matrix_postpc.loc[cell.id,2]]
		colors += [colors_dict[cell.library_number]]

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(comp_0, comp_1, comp_2, color = colors, s = 2, alpha = 1)

ax.set_xlabel('Component 0', fontsize = 12)
ax.set_ylabel('Component 1', fontsize = 12)
ax.set_zlabel('Component 2', fontsize = 12)
ax.set_xlim([-80,40])
ax.set_ylim([-30,50])
ax.set_zlim([-30,30])
ax.set_title('PCA of Kallisto counts', y= 1.1, fontsize = 14)
ax.legend(handles = [blue_line, red_line, green_line, black_line, cyan_line], loc = 2)
plt.savefig("plots/qc10_3dpca_color_by_library_number_1to5.pdf")

"""
Plot 3D pca components and label by library number (libraries 1-5)
"""

times = [0,75,150,300]
colors_dict = {6:'b',7:'r',8:'g',9:'k',10:'c'}

blue_line = mpl.lines.Line2D([],[],color = 'blue', label = 'Library 6')
red_line = mpl.lines.Line2D([],[],color = 'red', label = 'Library 7')
green_line = mpl.lines.Line2D([],[],color = 'green', label = 'Library 8')
black_line = mpl.lines.Line2D([],[],color = 'black', label = 'Library 9')
cyan_line = mpl.lines.Line2D([],[],color = 'cyan', label = 'Library 10')


comp_0 = []
comp_1 = []
comp_2 = []
colors = []

for cell in all_cells:
	if cell.library_number in [6,7,8,9,10]:
		comp_0 += [data_matrix_postpc.loc[cell.id,0]]
		comp_1 += [data_matrix_postpc.loc[cell.id,1]]
		comp_2 += [data_matrix_postpc.loc[cell.id,2]]
		colors += [colors_dict[cell.library_number]]

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(comp_0, comp_1, comp_2, color = colors, s = 2, alpha = 1)
ax.set_xlim([-80,40])
ax.set_ylim([-30,50])
ax.set_zlim([-30,30])
ax.set_xlabel('Component 0', fontsize = 12)
ax.set_ylabel('Component 1', fontsize = 12)
ax.set_zlabel('Component 2', fontsize = 12)

ax.set_title('PCA of Kallisto counts', y= 1.1, fontsize = 14)
ax.legend(handles = [blue_line, red_line, green_line, black_line, cyan_line], loc = 2)
plt.savefig("plots/qc10_3dpca_color_by_library_number_6to10.pdf")






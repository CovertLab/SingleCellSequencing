"""
Analysis!

Scatter plot pairs of genes
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
# import h5py
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
import scipy.stats as stats

from dba import dba
from dba import align_to

from rpy2.robjects.vectors import DataFrame as RDataFrame
from rpy2 import rinterface
from rpy2.robjects import conversion
from sklearn import manifold

@conversion.py2ro.register(pd.DataFrame)
def py2ro_pandasdataframe(obj):
    ri_dataf = conversion.py2ri(obj)
    # cast down to an R list (goes through a different code path
    # in the DataFrame constructor, avoiding `str(k)`) 
    ri_list = rinterface.SexpVector(ri_dataf)
    return RDataFrame(ri_list)

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = .1
# mpl.style.use('ggplot')

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')
scde = importr("scde")

# Load data sets in R
r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")
r("""load("/scratch/PI/mcovert/dvanva/sequencing/counts_data.RData")""")
r("o.fpm = scde.expression.magnitude(o.ifm, counts = counts_data_int)")
fpm = pd.read_csv('/scratch/PI/mcovert/dvanva/sequencing/fpms.csv', index_col = 0)
logfpm = np.log2(fpm + 1)

list_of_all_genes = list(r("rownames(counts_data_int)"))

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

time_points = [75, 150, 300]
cluster_list = [0,1,2]
color_dict = {'0':'g', '1':'r', '2':'b'}
# list_of_genes = ['Ccl3', 'Ccl4', 'Ccl5', 'Cxcl10', 'Il27','Saa3']
list_of_genes = ['Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl6', 'Ccl7', 'Ccl9', 'Ccl11', 'Ccl17', 'Ccl20', 'Ccl22', 'Ccl24', 'Ccl25', 'Ccl28']
list_of_genes += ['Cxcl1', 'Cxcl2', 'Cxcl3', 'Cxcl5', 'Cxcl9', 'Cxcl10', 'Cxcl11', 'Cxcl14', 'Cxcl15', 'Cxcl16', 'Cxcl17']
for time_point in time_points:
	all_cells = []
	cell_names = []
	longest_time = 0

	for cell in all_cells_total:
		if cell.time_point == time_point and cell.condition == 'Stim':
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]
			cell_names += [cell.id]

	"""
	Perform hierarchical clustering of the dynamics
	"""

	number_of_cells = len(all_cells)
	print number_of_cells
	distance_matrix_dynamics = np.zeros((number_of_cells, number_of_cells))

	dynamics_load = np.load('/home/dvanva/SingleCellSequencing/' + str(time_point)+'_dynamics_distance_matrix_kshape.npz')
	distance_matrix_dynamics = dynamics_load["distance_matrix"]

	Y_dynamics = sch.linkage(distance_matrix_dynamics, method = 'ward')
	ind_dynamics = sch.fcluster(Y_dynamics,0.5*np.amax(Y_dynamics[:,2]),'distance')
	
	ind_dynamics -= 1

	for j in xrange(number_of_cells):
		all_cells[j].clusterID = ind_dynamics[j]

		if time_point == 150:
			if ind_dynamics[j] == 2:
				all_cells[j].clusterID = 0
			if ind_dynamics[j] == 0:
				all_cells[j].clusterID = 1
			if ind_dynamics[j] == 1:
				all_cells[j].clusterID = 2

		for cell in all_cells_total:
			if cell.id == all_cells[j].id:
				cell.clusterID = all_cells[j].clusterID

for time_point in time_points:
	all_cells = []
	cell_names = []
	longest_time = 0

	plt.clf()
	fig, axes = plt.subplots(1,1, figsize = (6,6), squeeze = False)

	tsne = manifold.TSNE(n_components = 2, init = 'pca', random_state = 0)
	for cell in all_cells_total:
		if cell.time_point == time_point and cell.condition == 'Stim':
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]
			cell_names += [cell.id]

	reduced_fpm = logfpm.loc[list_of_genes,cell_names]
	trans_data = tsne.fit_transform(reduced_fpm.transpose()).T
	print trans_data.shape
	axes[0,0].scatter(trans_data[0], trans_data[1], s = 1, color = 'b')


	# for cluster in cluster_list:
	# 	print cluster
	# 	cell_names = []
	# 	for cell in all_cells:
	# 		if cell.clusterID == cluster:
	# 			cell_names += [cell.id]

	# 	reduced_fpm = logfpm.loc[:,cell_names]

	# 	counter1 = 0
	# 	for gene1 in list_of_genes:
	# 		counter2 = 0
	# 		for gene2 in list_of_genes:
	# 			gene1_array = reduced_fpm.loc[gene1,:]
	# 			gene2_array = reduced_fpm.loc[gene2,:]
	# 			axes[counter1,counter2].scatter(gene1_array, gene2_array, s = 3, color = color_dict[str(cluster)])
	# 			axes[counter1,counter2].set_xlabel(gene1)
	# 			axes[counter1,counter2].set_ylabel(gene2)

	# 			counter2 += 1
	# 		counter1 += 1

	fig.tight_layout()
	plt.savefig("plots/trial_31_tsne_" + str(time_point) + "min.pdf")










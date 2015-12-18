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
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2 import robjects as ro

pandas2ri.activate()
rpy2.robjects.numpy2ri.activate()
matplotlib.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
# matplotlib.style.use('ggplot')

# import pandas.DataFrame as PandasDataFrame
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

def extract_dataframe_from_R(dataframe_name):
	temp = pandas2ri.ri2py(r(dataframe_name))
	temp_rows = pandas2ri.ri2py(r("rownames(" + dataframe_name + ")"))
	temp_cols = np.float32(pandas2ri.ri2py(r("colnames(" + dataframe_name + ")")))

	df = pd.DataFrame(data = temp, columns = temp_cols, index = temp_rows)
	return df

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


"""
Load cells
"""
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_all_detected_genes_qc_w_jackpot.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Select and cluster 300 min time point
"""

time_point = 300

print "Analyzing " + str(time_point) + " minute time point"
all_cells = []
longest_time = 0
number_of_cells = 0

for cell in all_cells_total:
	if cell.time_point == time_point and cell.condition == 'Stim':
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
		all_cells += [cell]

all_cells_all_times = []
for t in [0,75,150,300]:
	for cell in all_cells_total:
		if cell.time_point == t and cell.condition == 'Stim':
			all_cells_all_times += [cell]
all_cells_total = all_cells_all_times

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

for cell in all_cells_total:
	if cell.id in cluster_dict:
		cell.clusterID = cluster_dict[cell.id]
	else:
		cell.clusterID = 0

"""
Send 300 min time point over to SCDE using R2py
"""

fact = ''
for cell in all_cells_total:
	fact += str(cell.clusterID)

fact_r = ro.FactorVector(fact)

r("""facts = """ + fact_r.r_repr())
r("""names(facts) = colnames(counts_data_int)""")

r("facts_3v1 = facts")
r("""facts_3v1[ facts == 0 | facts == 2] = NA""")
r("facts_3v1 = factor(facts_3v1)")
r("print(facts_3v1)")

r("facts_3v2 = facts")
r("""facts_3v2[ facts == 0 | facts == 1] = NA""")
r("facts_3v2 = factor(facts_3v2)")
r("print(facts_3v2)")

r("valid.cells = o.ifm$corr.a > 0")
r("table(valid.cells)")
r("o.ifm = o.ifm[valid.cells, ]")

r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")

r("ediff_3v1 = scde.expression.difference(o.ifm, counts_data_int, o.prior, groups  =  facts_3v1, n.randomizations  =  100, n.cores  =  1, verbose  =  1)")
r("""write.table(ediff_3v1[order(abs(ediff_3v1$Z), decreasing = TRUE), ], file = "/scratch/PI/mcovert/dvanva/sequencing/300min_scde_results_3v1.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)""")

r("ediff_3v2 = scde.expression.difference(o.ifm, counts_data_int, o.prior, groups  =  facts_3v2, n.randomizations  =  100, n.cores  =  1, verbose  =  1)")
r("""write.table(ediff_3v2[order(abs(ediff_3v2$Z), decreasing = TRUE), ], file = "/scratch/PI/mcovert/dvanva/sequencing/300min_scde_results_3v2.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)""")







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

# Load data sets in R
r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")
r("""load("/scratch/PI/mcovert/dvanva/sequencing/counts_data.RData")""")

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

# Determine which genes to look at
inflammatory_genes = ["Cxcl3", "Cxcl2", "Lif", "Ccl4", "Csf3", "Il1f9", "Ccl3", "Ccl5", "Tnf", "Il1a", "Il1b", "Tnfsf9", "Ccl20", "Il1f6", "Il27", "Il6"]
regulatory_genes = ["Nlrp3", "Nfkbiz", "Tnfaip2", "Nfkbia", "Tnfaip3", "Nfatc1"]
metabolic_genes = ["Hmox", "Prdx1", "Hdc", "Ptgs2", "Irg1"]
other_genes = ["Plaur", "Sqstm1", "Clec4e", "Sdc4", "Procr", "Slpi", "Plk2", "Saa3", "Slc7a11", "Cish", "Gp49a", "Hcar2", "Gpr84", "Malt1"]
inflammatory_genes = ["Ccl5"]

"""
Analyze all the time points
"""
cluster_list = {}
cluster_name_dict = {'0':{}, '75':{}, '150':{}, '300':{}}
times_to_analyze = [0, 300]
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
	if time_point != 0:
		dynamics_load = np.load('/home/dvanva/SingleCellSequencing/' + str(time_point)+'_dynamics_distance_matrix_kshape.npz')
		distance_matrix_dynamics = dynamics_load["distance_matrix"]

		Y_dynamics = sch.linkage(distance_matrix_dynamics, method = 'ward')
		ind_dynamics = sch.fcluster(Y_dynamics,0.5*np.amax(Y_dynamics[:,2]),'distance')

	if time_point == 0:
		cluster_list[str(time_point)] = np.arange(1,2)
	else:
		cluster_list[str(time_point)] = np.arange(np.amin(ind_dynamics), np.amax(ind_dynamics)+1)


	if time_point == 0:
		for j in xrange(number_of_cells):
			all_cells[j].clusterID = 1
	else:
		for j in xrange(number_of_cells):
			all_cells[j].clusterID = ind_dynamics[j]

	cluster_dict = {}

	for cell in all_cells:
		cluster_dict[cell.id] = str(cell.clusterID)

	for cluster in cluster_list[str(time_point)]:
		cluster_name_dict[str(time_point)][str(cluster)] = []
		for cell in all_cells:
			if cell.clusterID == cluster:
				cluster_name_dict[str(time_point)][str(cluster)] += [cell.id]

	"""
	Compute posterior FPM distribution for a given gene
	"""

	for gene in inflammatory_genes:
		gene_name = """'""" + gene + """'"""

		scde = importr("scde")

		r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
		r("""gene_counts = counts_data_int[c(""" + gene_name + ""","mt-Atp8"),]""")

		fpm_list = []
		jp_list = []
		ratio_list = []
		post_list = []
		for cluster in cluster_list[str(time_point)]:
			if time_point == 0:
				list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(time_point)][str(cluster)])
				r("list_of_cells = " + list_of_cells_r.r_repr())
				r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
				r("jp_0 = joint_posterior[" + gene_name + ",]")
				
				fpms = ro.r("colnames(joint_posterior)")
				fpms = np.float32(pandas2ri.ri2py(fpms))

				jp = ro.r("jp_0")
				jp = np.float32(pandas2ri.ri2py(jp))

				fpm_list += [fpms]
				jp_list += [jp]

			else:
				list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(time_point)][str(cluster)])
				r("list_of_cells = " + list_of_cells_r.r_repr())
				r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
				r("jp = joint_posterior[" + gene_name + ",]")
				
				fpms = ro.r("colnames(joint_posterior)")
				fpms = np.float32(pandas2ri.ri2py(fpms))

				jp = ro.r("jp")
				jp = np.float32(pandas2ri.ri2py(jp))

				fpm_list += [fpms]
				jp_list += [jp]

				r("ratio = scde:::calculate.ratio.posterior(t(jp), t(jp_0), o.prior, n.cores = 2)")

				ratios = 10 ** np.float32(pandas2ri.ri2py(r("colnames(ratio)")))
				post = np.float32(pandas2ri.ri2py(r("ratio")))
				ratio_list += [ratios]
				post_list += [post]

		"""
		Plot posteriors
		"""

		if time_point == 300:
			colors = ['g', 'r', 'b', 'k']

			plt.clf()

			max_jp = np.amax(jp_list[0])
			for j in xrange(len(fpm_list)):
				fpm = fpm_list[j]
				fpm_log2 = np.log2(fpm + 1e-50)
				jp = jp_list[j]

				print ratio_list[j].shape
				print post_list[j].shape
				max_jp = np.maximum(max_jp, np.amax(post_list[j]))
				plt.plot(ratio_list[j], np.squeeze(post_list[j]), color = colors[j], linewidth = 2, label = 'Cluster ' + str(j+1))
				plt.xlabel('log2(FPM)', fontsize = 16)
				plt.ylabel('Probability density', fontsize = 16)
				plt.title(gene + " " + str(time_point) + " minutes", fontsize = 16)
				plt.xlim([0,300])
				plt.xticks([0,300],  fontsize = 16)
				plt.ylim([0, 1.05*max_jp])
				plt.yticks([0, 1.05*max_jp],  fontsize = 16)

			plt.tight_layout()
			file_name = "trial_18_" + gene + "_" + str(time_point) + "min_ratio" + ".pdf"
			plt.savefig("plots/" + file_name)



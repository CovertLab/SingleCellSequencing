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
all_cell_file = 'all_cells_all_detected_genes_qc_w_jackpot.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

# Determine which genes to look at
inflammatory_genes = ["Cxcl3", "Cxcl2", "Lif", "Ccl4", "Csf3", "Il1f9", "Ccl3", "Ccl5", "Tnf", "Il1a", "Il1b", "Tnfsf9", "Ccl20", "Il1f6", "Il27", "Il6"]
regulatory_genes = ["Nlrp3", "Nfkbiz", "Tnfaip2", "Nfkbia", "Tnfaip3", "Nfatc1"]
metabolic_genes = ["Hmox", "Prdx1", "Hdc", "Ptgs2", "Irg1"]
other_genes = ["Plaur", "Sqstm1", "Clec4e", "Sdc4", "Procr", "Slpi", "Plk2", "Saa3", "Slc7a11", "Cish", "Gp49a", "Hcar2", "Gpr84", "Malt1"]
tnf_genes = ["Ptgs2", "Edn1", "Ereg", "Gadd45b" , "Ier2", "Cxcl1", "Cxcl2", "Il6", "Lif", "Btg2", "Zfp36", "Tnfaip3", "Ier3", "Sqstm1","Atf3","Ifrd1","Nfkbia","Nfkbiz","Irf1", "Junb","Slc25a25","Ifi47", "Icam1","Vcam1","Birc2", "Ccl2", "Ccl20", "Ccl7", "Cx3cl1","Cxcl10","Cxcl5", "Ifngr2", "Ripk2", "S1pr3","Tlr2", "Cebpd", "Nfkbie", "Relb", "Tnfaip2", "Rnd1", "Tmem140", "Stx11", "H2-D1", "Psmb10", "Tapbp", "C3", "Gbp2", "Ifit1","Casp4","Cdsn", "Mt2", "Sod3","Serpina3g","Ubd","Ccl5","Ccl9", "Lcn2","Il13ra1","Il1rl1","Naf1", "Enpp2","Mmp13", "Mmp3" "Serpinf1","Ppicap", "Ifitm3","Abca1", "Lhfpl2", "Stap2","Prnd", "U90926"]

list_of_genes = set(all_cells_total[0].fpkm.index)
tnf_genes = set(tnf_genes)
tnf_genes &= list_of_genes
tnf_genes = list(tnf_genes)
print len(tnf_genes)
print tnf_genes

"""
Analyze all the time points
"""

times_to_analyze = [75, 150, 300]
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

	"""
	Compute posterior FPM distribution for a given gene
	"""
	plt.clf()
	fig, axes = plt.subplots(len(tnf_genes)/4 + 1 ,4, figsize = (40,40))
	counter = 0

	for gene in tnf_genes:
		gene_name = """'""" + gene + """'"""

		scde = importr("scde")

		r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
		r("""gene_counts = counts_data_int[c(""" + gene_name + ""","mt-Atp8"),]""")

		fpm_list = []
		jp_list = []
		for cluster in cluster_list:
			list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(cluster)])
			r("list_of_cells = " + list_of_cells_r.r_repr())
			r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
			r("jp = joint_posterior[" + gene_name + ",]")
			
			fpms = ro.r("colnames(joint_posterior)")
			fpms = np.float32(pandas2ri.ri2py(fpms))

			jp = ro.r("jp")
			jp = np.float32(pandas2ri.ri2py(jp))

			fpm_list += [fpms]
			jp_list += [jp]


		"""
		Plot posteriors
		"""

		colors = ['g', 'r', 'b', 'k']

		max_jp = np.amax(jp_list[0])
		for j in xrange(len(fpm_list)):
			fpm = fpm_list[j]
			fpm_log2 = np.log2(fpm + 1e-50)
			jp = jp_list[j]

			max_jp = np.maximum(max_jp, np.amax(jp))
			axes.flatten()[counter].plot(fpm_log2, jp, color = colors[j], linewidth = 2, label = 'Cluster ' + str(j+1))
			axes.flatten()[counter].set_xlabel('log2(FPM)', fontsize = 16)
			axes.flatten()[counter].set_ylabel('Probability density', fontsize = 16)
			axes.flatten()[counter].set_title(gene + " " + str(time_point) + " minutes", fontsize = 16)
			axes.flatten()[counter].set_xlim([0,30])
			axes.flatten()[counter].set_xticks([0,10,20,30])
			axes.flatten()[counter].set_ylim([0, 1.05*max_jp])
			axes.flatten()[counter].set_yticks([0, 1.05*max_jp])
		counter += 1

	fig.tight_layout()
	file_name = "trial_19_all_genes_" + str(time_point) + "min" + ".pdf"
	plt.savefig("plots/" + str(time_point) + "min/tnf/" + file_name)



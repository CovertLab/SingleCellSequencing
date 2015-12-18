

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
import seaborn as sns
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

def zero_crossing(data, offset):
	new_data = data - offset
	zc = np.where(np.diff(np.signbit(new_data)))[0]
	return zc[0]

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

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

# Determine which genes to look at
# list_of_genes = ["Nfkbia", "Nfkbie", "Nfkbiz", "Nlrp3", "Tnfaip2", "Tnfaip3"]
# list_of_genes = ["Cxcl2", "Cxcl3", "Cxcl10", "Ccl3", "Ccl4", "Ccl5", "Ccl20", "Tnf", "Tnfsf9", "Il1a", "Il1b", "Il1f6", "Il6", "Il27", "Il1f9", "Lif",  "Csf3"]
# list_of_genes = ["Hmox1", "Prdx1", "Hdc", "Ptgs2", "Irg1"]
# list_of_genes = ["Plaur", "Sqstm1", "Clec4e", "Sdc4", "Procr", "Slpi", "Plk2", "Saa3", "Slc7a11", "Cish", "Gp49a", "Hcar2", "Gpr84", "Malt1"]
# list_of_genes = ['Dnmt3b', 'Tecta', 'Tm6sf2', 'Bricd5', 'Prdm12', 'Prdm13', 'Adora2a', 'Ccdc162', 'Gm5283', 'Gm11400', 'Olfr536', 'Gm13145', 'Gm13333', 'Zfp661', 'Angptl3', 'Sipa1l3', 'Scn1a', 'Sprr2d', 'Il17rc', 'Zglp1', 'Akr1cl', 'Map1a', 'Trim14', 'Adgrg6', 'Gm13991', 'Dhrs11', 'Gm21834', 'Iqca', 'Gm2007', 'Slc39a8', 'Gng7', 'AL663030.1', 'Nphp4', 'Nod1', 'Emc9', 'Akr1b7', 'Il33', 'Mmp14', 'Zfyve1', 'Cetn4', '2610305D13Rik', 'Mettl25', 'Ric8b', 'Mterf2', 'Zfp850', 'Clec4a4', 'Saa3', 'Hist1h4n', 'Gm11007', 'Cntrob', 'Atp7b', 'Mtl5', '1700061G19Rik', 'Coro2b', '1700030J22Rik', 'Gm8898', 'Tmem86b', 'Car9', 'Gm5157', 'Gm15539', 'Arhgef18', 'Slc13a3', 'Dclk1', 'Ager', 'Actr3b', 'Zfp41', 'Fzd8', '4930524J08Rik', 'Zic5', 'Trem1', 'Ppp1r32', 'Stk36', 'Gnao1', 'Tmem239', 'Polm', 'Fgf21', 'Gprasp2', 'Tesk1', 'Athl1', 'Kptn']
# list_of_genes = ['Mmp3', 'Ccl5', 'Gpr137c', 'Efna5', 'Tiam1', 'D2hgdh', 'Nod2', 'Gm14440', 'Pla2r1', 'Serpinb9g', 'Hic2', 'Cdkl4', 'Slc18b1', 'H2-M2', 'Klhdc1', 'Iqcb1', 'Sh3bp2', 'Ifit3', 'Cmpk2', 'Adamts10', 'Sirt5', 'Plekhg2', 'Cxcl10', 'Gm13051', 'Tppp3', 'Krt24', 'Lamb3', 'Serpind1', 'Pars2', 'Spopl', 'Rsad2', 'Tnfsf4', 'Gm12728', 'Siglece', '4930432K21Rik', 'Vmn1r32', 'Fbxw10', 'Ngb', 'Bdkrb1', 'B3galt2']
list_of_genes = ["Ccl3"]

"""
Analyze all the time points
"""

times_to_analyze = [0, 75, 150, 300]
cluster_list = {}
cluster_name_dict = {'0':{}, '75':{}, '150':{}, '300':{}}

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

plt.clf()
# fig, axes = plt.subplots(len(list_of_genes)/4,4, figsize = (60,60))
counter = 0

for gene in list_of_genes:
	print gene
	cluster_1_mean = []
	cluster_2_mean = []
	cluster_3_mean = []

	cluster_1_low = []
	cluster_2_low = []
	cluster_3_low = []

	cluster_1_high = []
	cluster_2_high = []
	cluster_3_high = []

	for time_point in times_to_analyze:

		gene_name = """'""" + gene + """'"""

		r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
		r("""gene_counts = counts_data_int[c(""" + gene_name + ""","mt-Atp8"),]""")

		ratio_list = []
		post_list = []


		for cluster in cluster_list[str(time_point)]:
			if time_point == 0:
				list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(time_point)][str(cluster)])
				r("list_of_cells = " + list_of_cells_r.r_repr())
				r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
				r("prior = o.prior")
				r("jp_0 = joint_posterior[" + gene_name + ",]")
				r("jp_0 = t(jp_0)")
				r("ratio = scde:::calculate.ratio.posterior(jp_0, jp_0, prior = o.prior, n.cores = 2)")

				ratios = 10 ** np.float32(pandas2ri.ri2py(r("colnames(ratio)")))
				post = np.float32(pandas2ri.ri2py(r("ratio")))

				ratio_list += [ratios]
				post_list += [post]

			else:
				list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(time_point)][str(cluster)])
				r("list_of_cells = " + list_of_cells_r.r_repr())
				r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
				r("jp = joint_posterior[" + gene_name + ",]")
				r("jp = t(jp)")
				r("ratio = scde:::calculate.ratio.posterior(jp, jp_0, o.prior, n.cores = 2)")

				ratios = 10 ** np.float32(pandas2ri.ri2py(r("colnames(ratio)")))
				post = np.float32(pandas2ri.ri2py(r("ratio")))

				ratio_list += [ratios]
				post_list += [post]

		# Give the clusters the proper order
		if time_point == 150:
			ratio_list_new = []
			post_list_new = []

			ratio_list_new += [ratio_list[2]]
			ratio_list_new += [ratio_list[0]]
			ratio_list_new += [ratio_list[1]]
			post_list_new += [post_list[2]]
			post_list_new += [post_list[0]]
			post_list_new += [post_list[1]]
			
			ratio_list = ratio_list_new
			post_list = post_list_new

		if time_point == 0:
			ratio = ratio_list[0]
			post = post_list[0]

			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]

			cluster_1_mean += ratios
			cluster_2_mean += ratios
			cluster_3_mean += ratios

			cluster_1_low += err_low
			cluster_2_low += err_low
			cluster_3_low += err_low

			cluster_1_high += err_high
			cluster_2_high += err_high
			cluster_3_high += err_high

		if time_point == 75:
			ratio= ratio_list[0]
			post = post_list[0]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_1_mean += ratios
			cluster_3_mean += ratios
			cluster_1_low += err_low
			cluster_3_low += err_low
			cluster_1_high += err_high
			cluster_3_high += err_high

			ratio = ratio_list[1]
			post = post_list[1]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_2_mean += ratios
			cluster_2_low += err_low
			cluster_2_high += err_high

		if time_point == 150:
			ratio = ratio_list[0]
			post = post_list[0]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_1_mean += ratios
			cluster_1_low += err_low
			cluster_1_high += err_high

			ratio = ratio_list[1]
			post = post_list[1]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_2_mean += ratios
			cluster_2_low += err_low
			cluster_2_high += err_high

			ratio = ratio_list[2]
			post = post_list[2]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_3_mean += ratios
			cluster_3_low += err_low
			cluster_3_high += err_high

		if time_point == 300:
			ratio = ratio_list[0]
			post = post_list[0]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_1_mean += ratios
			cluster_1_low += err_low
			cluster_1_high += err_high

			ratio = ratio_list[1]
			post = post_list[1]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_2_mean += ratios
			cluster_2_low += err_low
			cluster_2_high += err_high

			ratio = ratio_list[2]
			post = post_list[2]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_3_mean += ratios
			cluster_3_low += err_low
			cluster_3_high += err_high

	"""
	Plot posteriors
	"""

	colors = ['g', 'r', 'b']

	cluster_1_low = np.array(cluster_1_low)
	cluster_2_low = np.array(cluster_2_low)
	cluster_3_low = np.array(cluster_3_low)

	cluster_1_high = np.array(cluster_1_high)
	cluster_2_high = np.array(cluster_2_high)
	cluster_3_high = np.array(cluster_3_high)

	cluster_1_mean = np.array(cluster_1_mean)
	cluster_2_mean = np.array(cluster_2_mean)
	cluster_3_mean = np.array(cluster_3_mean)

	cluster_means = np.zeros((3,4), dtype = 'float32')
	cluster_means[0,:] = cluster_2_mean
	cluster_means[1,:] = cluster_3_mean
	cluster_means[2,:] = cluster_1_mean

	# print [np.abs(cluster_3_low - cluster_3_mean), np.abs(cluster_3_high - cluster_3_mean)]
	# max_val = np.amax(np.array([np.amax(cluster_1_high), np.amax(cluster_2_high), np.amax(cluster_3_high)]))
	# axes.flatten()[counter].errorbar(times_to_analyze, cluster_1_mean, yerr = [np.abs(cluster_1_low - cluster_1_mean), np.abs(cluster_1_high - cluster_1_mean)], fmt = '-o', color = colors[0], ecolor = colors[0], linewidth = 2, label = 'Cluster 1')
	# axes.flatten()[counter].errorbar(times_to_analyze, cluster_2_mean, yerr = [np.abs(cluster_2_low - cluster_2_mean), np.abs(cluster_2_high - cluster_2_mean)], fmt = '-o', color = colors[1], ecolor = colors[1], linewidth = 2, label = 'Cluster 2')
	# axes.flatten()[counter].errorbar(times_to_analyze, cluster_3_mean, yerr = [np.abs(cluster_3_low - cluster_3_mean), np.abs(cluster_3_high - cluster_3_mean)], fmt = '-o', color = colors[2], ecolor = colors[2], linewidth = 2, label = 'Cluster 3')

	# axes.flatten()[counter].set_xlabel('Time (minutes)', fontsize = 16)
	# axes.flatten()[counter].set_ylabel('Fold change', fontsize = 16)
	# axes.flatten()[counter].set_title(gene, fontsize = 16)
	# axes.flatten()[counter].set_ylim([0,np.ceil(1.05*max_val)])
	# axes.flatten()[counter].set_yticks([0,np.ceil(1.05*max_val)])
	# axes.flatten()[counter].set_xlim([0, 1.05*np.amax(times_to_analyze)])
	# axes.flatten()[counter].set_xticks(times_to_analyze)
	ax = sns.heatmap(cluster_means, linewidths = 3, cmap = sns.light_palette("muted purple", input = "xkcd", as_cmap = True), vmin = 1, vmax = np.amax(cluster_means.flatten()), xticklabels = False, yticklabels = [2,3,1])
	ax.set_title(gene)
	counter += 1

# fig.tight_layout()
file_name = "trial_25_heatmap.pdf"
plt.savefig("plots/" + file_name)



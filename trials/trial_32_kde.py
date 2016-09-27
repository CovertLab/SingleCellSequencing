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
from scipy import interpolate
from dba import dba
from dba import align_to

from rpy2.robjects.vectors import DataFrame as RDataFrame
from rpy2 import rinterface
from rpy2.robjects import conversion

def zero_crossing(data, offset):
	new_data = data - offset
	zc = np.where(np.diff(np.signbit(new_data)))[0]
	return zc[0]

def randomvariate(pdf,n=100, xmin=0, xmax = 1):
	x = np.linspace(xmin,xmax,1000)
	y = pdf(x)
	pmin = 0.
	pmax = y.max()

	# Counters
	naccept = 0 
	ntrial = 0

	# Keep generating numbers until we achieve the desired n
	ran = []
	while naccept < n:
		x = np.random.uniform(xmin,xmax)
		y = np.random.uniform(pmin,pmax)

		if y < pdf(x):
			ran.append(x)
			naccept += 1
			ntrial += 1

	ran = np.asarray(ran)
	return ran

def get_individual_prior(cell_name, gene):
	list_of_cells_r = ro.vectors.StrVector([cell_name])
	gene_name = gene_name = """'""" + gene + """'"""

	r("list_of_cells = " + list_of_cells_r.r_repr())
	r("""gene_counts = counts_data_int[c(""" + gene_name + ""","mt-Atp8"),]""")
	r("""ind_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
	r("ip = ind_posterior[" + gene_name + ",]")
	
	fpms = ro.r("colnames(ind_posterior)")
	fpms = np.float32(pandas2ri.ri2py(fpms))

	ip = ro.r("ip")
	ip = np.float32(pandas2ri.ri2py(ip))

	return fpms, ip

def get_joint_posterior(cell_names, gene):
	list_of_cells_r = ro.vectors.StrVector(cell_names)
	gene_name = gene_name = """'""" + gene + """'"""

	r("list_of_cells = " + list_of_cells_r.r_repr())
	r("""gene_counts = counts_data_int[c(""" + gene_name + ""","mt-Atp8"),]""")
	r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
	r("jp = joint_posterior[" + gene_name + ",]")
	
	fpms = ro.r("colnames(joint_posterior)")
	fpms = np.float32(pandas2ri.ri2py(fpms))

	jp = ro.r("jp")
	jp = np.float32(pandas2ri.ri2py(jp))

	return fpms, jp

def interp_prior(cell_name = None, gene_name = None):
	fpms, ip = get_individual_prior(cell_name, gene_name)
	ip_rescale = ip/ip.max()
	loc = zero_crossing(ip_rescale, 0.005)

	limits = (np.log2(fpms[loc.min()]+1), np.log2(fpms[loc.max()]+1))

	return interpolate.interp1d(np.log2(fpms+1), ip), limits

def sample_prior(cell_name = None, gene_name = None):
	f, limits = interp_prior(cell_name, gene_name)
	ran = randomvariate(f, xmin = limits[0], xmax = limits[1])
	return ran

def computeMI(P):

	P = P/P.sum()
	px = P.sum(axis = 1)/P.sum(axis = 1).sum()
	py = P.sum(axis = 0)/P.sum(axis = 0).sum()


	px_tile = np.tile(px, (P.shape[1], 1))
	py_tile = np.tile(py, (P.shape[1],1)).T

	Px = P/px_tile
	Px /= Px.max(axis = 0)

	temp = Px * np.log2(P/(px_tile*py_tile))
	return temp.sum()


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

r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
r("o.fpm = scde.expression.magnitude(o.ifm, counts = counts_data_int)")


fpm = pd.read_csv('/scratch/PI/mcovert/dvanva/sequencing/fpms.csv', index_col = 0)
logfpm = np.log2(fpm + 1)

list_of_all_genes = list(r("rownames(counts_data_int)"))

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

time_points = [75,150,300]
color_dict = {'0':'g', '1':'r', '2':'b'}
# list_of_genes = ['Ccl3', 'Ccl4', 'Ccl5', 'Cxcl10', 'Il27','Saa3']
# list_of_genes = ['Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl6', 'Ccl7', 'Ccl9', 'Ccl11', 'Ccl17', 'Ccl20', 'Ccl22', 'Ccl24', 'Ccl25', 'Ccl28']
# list_of_genes += ['Cxcl1', 'Cxcl2', 'Cxcl3', 'Cxcl5', 'Cxcl9', 'Cxcl10', 'Cxcl11', 'Cxcl14', 'Cxcl15', 'Cxcl16', 'Cxcl17']
# list_of_genes = ['Ccl3', 'Ccl4', 'Ccl5', 'Cxcl10']
# list_of_genes = ['Nupr1', 'Cxcl10', 'Ccl5', 'Il27', 'Lilrb4', 'Oasl1', 'Icam1', 'Nfkbie', 'Irg1', 'Marcksl1', 'Slfn2', 'Gp49a', 'Spp1', 'Ccl3', 'Ccl4', 'Plaur', 'Srgn', 'Pim1', 'Traf1', 'Btg2', 'Ier3', 'Phlda1', 'Kdm6b', 'Nfkbia', 'Nfkbiz', 'Tfec', 'Malt1', 'Nlrp3', 'Sdc4', 'Cd83', 'Tnfaip3', 'Il1a', 'Il1b', 'Sat1', 'Cxcl3', 'Saa3', 'Sod2', 'Csf3', 'Cxcl2', 'Tnf', 'Il1f9', 'Tnfsf9']
# list_of_genes = ['Ccl3', 'Ccl4', 'Ccl5', 'Cxcl10', 'Saa3', 'Cxcl2', 'Nupr1', 'Irg1', 'Tnf', 'Nfkbia', 'Nfkbie', 'Nfkbiz', 'Tnfaip3','Dynamics']
list_of_genes = ['Ccl3', 'Ccl4', 'Ccl5', 'Cxcl10', 'Saa3', 'Cxcl2', 'Tnf', 'Nfkbia', 'Tnfaip3','Dynamics']

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
	fig, axes = plt.subplots(len(list_of_genes),len(list_of_genes), figsize = (4*len(list_of_genes),4*len(list_of_genes)))
	corr_matrix = np.zeros((len(list_of_genes), len(list_of_genes)))
	for cell in all_cells_total:
		if cell.time_point == time_point and cell.condition == 'Stim':
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]
			cell_names += [cell.id]

	all_cells_new = []
	if time_point == 75:
		cluster_list = [0,1]
	else:
		cluster_list = [0,1,2]
	for cluster in cluster_list:
		print cluster
		cell_names = []
		for cell in all_cells:
			if cell.clusterID == cluster:
				cell_names += [cell.id]
				all_cells_new += [cell]
	all_cells = all_cells_new

	for cell in all_cells:
		logfpm.loc['Dynamics', cell.id] = np.sum(cell.NFkB_dynamics)
		logfpm.loc['ClusterID', cell.id] = cell.clusterID
	reduced_fpm = logfpm.loc[list_of_genes,cell_names]


	# print reduced_fpm
	# print reduced_fpm.loc['ClusterID',:]
	counter1 = 0
	for gene1 in list_of_genes:
		counter2 = 0
		for gene2 in list_of_genes:
			print gene1, gene2
			# if counter1 > counter2:

			if gene1 == gene2:

				Z = np.rot90(np.eye(100))+1e-8
				corr_matrix[counter1, counter2]  = computeMI(Z)

				axes[counter1, counter2].imshow(Z, cmap = plt.get_cmap('coolwarm'), interpolation = 'none')

			if gene1 != gene2:
				# gene1_list = []
				# gene2_list = []
				# for cell_name in cell_names:


				# 	gene1_list += list(sample_prior(cell_name = cell_name, gene_name = gene1))
				# 	gene2_list += list(sample_prior(cell_name = cell_name, gene_name = gene2))

				# gene1_array = np.array(gene1_list)
				# gene2_array = np.array(gene2_list)

				# Non sampled implementation
				gene1_array = reduced_fpm.loc[gene1,:]
				gene2_array = reduced_fpm.loc[gene2,:]

				if gene1_array.sum() > 0 and gene2_array.sum() > 0:

					xmin = gene1_array.min()
					xmax = gene1_array.max()

					ymin = gene2_array.min()
					ymax = gene2_array.max()

					X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
					positions = np.vstack([X.ravel(), Y.ravel()])
					values = np.vstack([gene1_array, gene2_array])
					kernel = scipy.stats.gaussian_kde(values)
					Z = np.reshape(kernel(positions).T, X.shape)
					Z = np.rot90(Z)

					px = Z.sum(axis = 1)/Z.sum(axis = 1).sum()
					px_tile = np.tile(px, (Z.shape[1], 1))
					Zx = Z/px_tile
					Zx /= Zx.max(axis = 0)

					axes[counter1, counter2].imshow(Zx, cmap = plt.get_cmap('coolwarm'), interpolation = 'none', extent = [xmin, xmax, ymin, ymax])
					# axes[counter1,counter2].scatter(gene1_array, gene2_array, s = 3, color = 'b')
					axes[counter1,counter2].set_xlabel(gene1, fontsize = 14)
					axes[counter1,counter2].set_ylabel(gene2, fontsize = 14)

					# corr_matrix[counter1, counter2]  = computeMI(Z)
			counter2 += 1
		counter1 += 1

	print corr_matrix
	fig.tight_layout()
	plt.savefig("plots/trial_32_gene_gene_scatter_" + str(time_point) + "min.pdf")

	# Y = sch.linkage(corr_matrix, method = 'ward')
	# Z = sch.dendrogram(Y, orientation = 'right', color_threshold = 0.5*np.amax(Y[:,2]))
	# index = Z['leaves']
	# corr1 = corr_matrix[index, :]
	# corr_ordered = corr1[:,index]
	# list_of_genes_ordered = [list_of_genes[i] for i in index]

	# fig = plt.figure()
	# ax = fig.add_subplot(111)

	# ax.matshow(corr_matrix, aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('coolwarm'), interpolation = 'none') #, vmin = -1, vmax = 1)
	# ticks = np.arange(0,corr_matrix.shape[0],1)
	# ax.xaxis.set_ticks_position('bottom')
	# ax.set_xticks(ticks)
	# ax.set_yticks(ticks)
	# ax.set_xticklabels(list_of_genes, fontsize = 5, rotation = 270)
	# ax.set_yticklabels(list_of_genes, fontsize = 5)
	# ax.tick_params(width = 0, length = 0)

	# plt.savefig("plots/trial_32_MI_heatmap" + str(time_point) + "min.pdf")












	# cluster_dict = {}
	# cluster_name_dict = {}

	# for cell in all_cells:
	#   cluster_dict[cell.id] = str(cell.clusterID)

	# for cluster in cluster_list:
	#   cluster_name_dict[str(cluster)] = []
	#   for cell in all_cells:
	#       if cell.clusterID == cluster:
	#           cluster_name_dict[str(cluster)] += [cell.id]
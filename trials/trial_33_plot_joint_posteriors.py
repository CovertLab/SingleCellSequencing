"""
Analysis!

Plot joint posteriors for different clusters
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

def get_error_bar(fpms, jp, confidence_level = 0.16):
	cumsum = np.cumsum(jp)
	err_low = fpms[zero_crossing(cumsum, confidence_level)]
	err_high = fpms[zero_crossing(cumsum, 1-confidence_level)]
	return err_low, err_high

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = .1

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')
scde = importr("scde")

# Load data sets in R
r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")
r("""load("/scratch/PI/mcovert/dvanva/sequencing/counts_data.RData")""")

r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
r("o.fpm = scde.expression.magnitude(o.ifm, counts = counts_data_int)")

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))
time_points = [300]

# list_of_genes = ["Saa3", "Il1f9", "Ccl3", "Ccl5", "Prdx1"]

list_of_genes = ['Wrb', 'Cxcl3', 'Tnfaip3', 'Tnfaip2', 'Blvrb', 'Nfatc1', 'Nfkbiz', 'Phf11d', 'Fas', 'Gpx1', 'Ccnd1', 'Rnf213', 'Marcksl1', 'Abcg1', 'Gm8818', 'Nlrp3', 'Rrm2', 'Mbnl1', 'Ptpn14', 'Odc1', 'Ptgs2', 'Ddit3', 'Cxcl10', 'Ly86', 'Ier3', 'Cd14', 'Rel', 'Prkg2', 'Afp', 'Btg2', 'Gm18445', 'Sdc4', 'Tnfsf9', 'Prr14l', 'Il27', 'Tk1', 'Angpt2', 'Tmem171', 'Ccl3', 'Ccl4', 'Sqstm1', 'Cd83', 'Slc7a11', 'Oasl1', 'Hsp90aa1', 'Pde4b', 'Rasgrp3', 'Calcrl', 'Egr1', 'Stx11', 'Colec12', 'Gmnn', 'Gpr84', 'Cxcl2', 'Sod2', 'Mt2', 'Serpinb8', 'Srxn1', 'Phlda1', 'Bcl2a1d', 'Traf1', 'Pim1', 'Il1f9', 'Prdx1', 'Procr', 'Hmgcs1', 'AI607873', 'Bcl2l11', 'Irg1', 'Saa3', 'Tnf', 'Hdc', 'Atf3', 'Errfi1', 'Lif', 'Sat1', 'Plaur', 'Hmox1', 'Id1', 'Mef2c', 'Icam1', 'Slfn2', 'Map1b', 'Lilrb4', 'Clec4e', 'Nfkbia', 'Csf3', 'Akr1b8', 'Emp1', 'Srgn', 'Zc3h12c', 'Slpi', 'Rasgef1b', 'Plk3', 'Plk2', 'Rassf4', 'Stap1', 'Kdm6b', 'Il1b', 'Gp49a', 'Malt1', 'Nabp1', 'Kif14', 'Rab31', 'Ppp1r15a', 'Mtss1', 'Ccl5']

cluster_name_dict = {'0':{}, '75':{}, '150':{}, '300':{}}
colors = ['g', 'r', 'b', 'k']
plt.clf()
fig, axes = plt.subplots(1,1, figsize = (10,10), squeeze = False)

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

	if time_point == 0:
		cluster_list = np.array([0])
	else:
		cluster_list = np.array([0,1,2])

	for cluster in cluster_list:
		cluster_name_dict[str(time_point)][str(cluster)] = []
		for cell in all_cells:
			if cell.clusterID == cluster:
				cluster_name_dict[str(time_point)][str(cluster)] += [cell.id]
	
	cluster_fpms = {"0":[], "1":[], "2":[],"all":[]}
	cluster_err_high = {"0":[], "1":[], "2":[],"all":[]}
	cluster_err_low = {"0":[], "1":[], "2":[],"all":[]}


	for gene in list_of_genes:
		max_jp = 0
		for cluster in cluster_list:
			colors = ['g', 'r', 'b', 'k']
			cell_names = cluster_name_dict[str(time_point)][str(cluster)]
			fpm, jp = get_joint_posterior(cell_names, gene)
			max_jp = np.maximum(max_jp,np.amax(jp))
			fpm_log2 = np.log2(fpm + 1e-50)
			cluster_fpms[str(cluster)] +=  [fpm_log2[np.argmax(jp)]]
			err_low, err_high = get_error_bar(fpm, jp)
			cluster_err_high[str(cluster)] += [np.log2(err_high + 1e-50)]
			cluster_err_low[str(cluster)] += [np.log2(err_low + 1e-50)]

		cell_names = []
		for cluster in cluster_list:
			cell_names += cluster_name_dict[str(time_point)][str(cluster)]
		fpm, jp = get_joint_posterior(cell_names, gene)
		cluster_fpms["all"] +=  [fpm_log2[np.argmax(jp)]]
		err_low, err_high = get_error_bar(fpm, jp)
		cluster_err_high["all"] += [np.log2(err_high + 1e-50)]
		cluster_err_low["all"] += [np.log2(err_low + 1e-50)]
		max_jp = np.maximum(max_jp,np.amax(jp))

	for cluster in cluster_list:
		axes[0,0].scatter(cluster_fpms["all"], cluster_fpms[str(cluster)], s = 4, color = colors[cluster])


	x = np.arange(0,20,.05)
	y = np.arange(0,20,.05)

	axes[0,0].plot(x,y,color = 'k')
	axes[0,0].plot(x,y+1,color = 'k')
	axes[0,0].plot(x+1,y,color = 'k')

	axes[0,0].set_xlim([0, 20])
	axes[0,0].set_ylim([0, 20])

	for counter in xrange(len(list_of_genes)):
		for cluster in cluster_list:
			if list_of_genes[counter] == 'Ccl5' or list_of_genes[counter] == 'Cxcl10' or list_of_genes[counter] == 'Saa3':
				axes[0,0].annotate(list_of_genes[counter], xy = (cluster_fpms["all"][counter], cluster_fpms[str(cluster)][counter]), textcoords = 'data', fontsize = 3, color = 'r')
			else:
				axes[0,0].annotate(list_of_genes[counter], xy = (cluster_fpms["all"][counter], cluster_fpms[str(cluster)][counter]), textcoords = 'data', fontsize = 3, color = 'k')

	axes.flatten()[0].set_xlabel('log2(FPM+1) (no clustering)', fontsize = 16)
	axes.flatten()[0].set_ylabel('log2(FPM+1) (clustering with dynamics)', fontsize = 16)
	fig.tight_layout()
	file_name = "trial_33_" + str(time_point) + "min" + ".pdf"
	plt.savefig("plots/" + file_name)


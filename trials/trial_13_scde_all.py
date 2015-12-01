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

"""
Load cells
"""
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_all_detected_genes_qc_w_jackpot.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Select and cluster 75 min time point
"""

times = [0,75,150,300]

R = rpy2.robjects.r
DTW = importr('dtw')

all_cells = []
longest_time = 0
number_of_cells = 0
for t in times:
	for cell in all_cells_total:
		if cell.time_point == t and cell.condition == 'Stim':
			number_of_cells += 1
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]

print len(all_cells)


for j in xrange(number_of_cells):
	if cell.time_point == 0:
		all_cells[j].clusterID = str(0)
	if cell.time_point == 75:
		all_cells[j].clusterID = str(1)
	if cell.time_point == 150:
		all_cells[j].clusterID = str(2)
	if cell.time_point == 300:
		all_cells[j].clusterID = str(3)

list_of_cell_names = []
for cell in all_cells:
	list_of_cell_names += [cell.id]

print list_of_cell_names

list_of_genes = all_cells[0].tpm.index.get_values()
number_of_genes = len(list_of_genes)

"""
Send cells over to SCDE using Rpy2
"""

counts_data = pd.DataFrame(np.zeros((number_of_genes, number_of_cells), dtype = 'int16'), index = list_of_genes, columns = list_of_cell_names)
for cell in all_cells:
	total_transcript_counts = cell.total_estimated_counts - cell.spikeins.est_counts.sum()
	est_counts =  cell.fpkm.transpose() *cell.mean_eff_length.transpose() * total_transcript_counts/1e9  #cell.est_counts.transpose() #*total_transcript_counts/1e3
	counts_data.loc[:,cell.id] = np.floor(est_counts)

factors = pd.DataFrame(np.zeros((1,number_of_cells), dtype = 'int16'), columns = list_of_cell_names)
fact = ''
for cell in all_cells:
	fact += str(cell.clusterID)

scde = importr("scde")

counts_data_r = conversion.py2ro(counts_data)
fact_r = ro.FactorVector(fact)

counts_data_r.colnames = list_of_cell_names

r("""counts_data_int = apply(""" + counts_data_r.r_repr() + """,2, function(x) {storage.mode(x) = 'integer';x}) """)

r("""facts = """ + fact_r.r_repr())
r("""names(facts) = colnames(counts_data_int)""")

r("o.ifm = scde.error.models(counts = counts_data_int, groups = facts, n.cores = 2, linear.fit=F, local.theta.fit=F, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)")
r("""save(o.ifm, file = "/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")

r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")

r("valid.cells = o.ifm$corr.a > 0")
r("table(valid.cells)")
r("o.ifm = o.ifm[valid.cells, ]")

r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
r("o.fpm = scde.expression.magnitude(o.ifm, counts = counts_data_int)")
r("o.fail.curves = scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))")

"""
Compute differential expression for all time points compared to the zero time point
"""

r("""facts_75 = facts""")
r("""facts_75[ facts == "150" | facts == "300"] = NA""")

r("""facts_150 = facts""")
r("""facts_150[ facts == "75" | facts == "300"] = NA""")

r("""facts_300 = facts""")
r("""facts_300[ facts == "75" | facts == "150"] = NA""")

r("ediff = scde.expression.difference(o.ifm, counts_data_int, o.prior, groups  =  facts_75, n.randomizations  =  100, n.cores  =  1, verbose  =  1)")
r("""write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "/scratch/PI/mcovert/dvanva/sequencing/75min_scde_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)""")

r("ediff = scde.expression.difference(o.ifm, counts_data_int, o.prior, groups  =  facts_150, n.randomizations  =  100, n.cores  =  1, verbose  =  1)")
r("""write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "/scratch/PI/mcovert/dvanva/sequencing/150min_scde_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)""")

r("ediff = scde.expression.difference(o.ifm, counts_data_int, o.prior, groups  =  facts_300, n.randomizations  =  100, n.cores  =  1, verbose  =  1)")
r("""write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "/scratch/PI/mcovert/dvanva/sequencing/300min_scde_results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)""")

"""
Plot posteriors for different clusters
"""

# r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/75_min_tnf_diff.pdf")""")
# r("""scde.test.gene.expression.difference("Tnf", models = o.ifm, counts = counts_data_int, prior = o.prior)""")
# r("dev.off()")

# r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/75_min_cxcl2_diff.pdf")""")
# r("""scde.test.gene.expression.difference("Cxcl2", models = o.ifm, counts = counts_data_int, prior = o.prior)""")
# r("dev.off()")

# r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/75_min_ccl22_diff.pdf")""")
# r("""scde.test.gene.expression.difference("Ccl22", models = o.ifm, counts = counts_data_int, prior = o.prior)""")
# r("dev.off()")

# r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/75_min_ccl4_diff.pdf")""")
# r("""scde.test.gene.expression.difference("Ccl4", models = o.ifm, counts = counts_data_int, prior = o.prior)""")
# r("dev.off()")

# r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/75min_cxcl10_diff.pdf")""")
# r("""scde.test.gene.expression.difference("Cxcl10", models = o.ifm, counts = counts_data_int, prior = o.prior)""")
# r("dev.off()")

"""
Plot failure curves
"""

r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/75min_fail_curves.pdf")""")
r("par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)")
r("""plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")""")
r("""invisible(apply(o.fail.curves[, grep("0",colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))""")
r("""invisible(apply(o.fail.curves[, grep("1", colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))""")
r("dev.off()")

r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/150min_fail_curves.pdf")""")
r("par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)")
r("""plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")""")
r("""invisible(apply(o.fail.curves[, grep("0",colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))""")
r("""invisible(apply(o.fail.curves[, grep("2", colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))""")
r("dev.off()")

r("""pdf("/home/dvanva/SingleCellSequencing/trials/plots/300min_fail_curves.pdf")""")
r("par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)")
r("""plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")""")
r("""invisible(apply(o.fail.curves[, grep("0",colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))""")
r("""invisible(apply(o.fail.curves[, grep("3", colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))""")
r("dev.off()")



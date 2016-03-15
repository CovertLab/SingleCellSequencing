"""
smFISH data analysis
Perform analysis for single cell analysis data
"""

"""
Import python packages
"""

import HTSeq 
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
import scipy.cluster.hierarchy as sch
from seq_functions import smFISH_cell, cleanAxis
import rpy2
from rpy2.robjects.packages import importr
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import sem
from scipy.ndimage.filters import median_filter as med_filt
import matplotlib.ticker as mtick
from scipy import stats
import scikits.bootstrap as bootstrap
import seaborn as sns

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
# mpl.style.use('ggplot')

rpy2.robjects.numpy2ri.activate()

"""
Load cluster averages
"""
direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name = os.path.join(direc,"300_cluster_avg_kshape_smFISH.npz")
file_load = np.load(file_name)
cluster_dynamics_avg = file_load["cluster_dynamics_avg"]

"""
Load good cells
"""
direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name_save = os.path.join(direc, "good_cells_300min.pkl")
good_cells = pickle.load(open(file_name_save))

numbers = np.array([0,0,0])
for cell in good_cells:
	if cell.clusterID == 0:
		numbers[0] += 1
	if cell.clusterID == 1:
		numbers[1] += 1
	if cell.clusterID == 2:
		numbers[2] += 1

print numbers


"""
Plot individual traces
"""

# # plt.plot(cluster_dynamics_avg[2,:])
# fig, axes = plt.subplots(3,10, figsize = (20,10))
# print axes.shape
# for cluster in [0,1,2]:
# 	counter = 0

# 	for cell in good_cells:
# 		if cell.clusterID == cluster and counter < 10:	
# 			print cluster, counter	
# 			axes[cluster, counter].plot(cell.norm_med)
# 			axes[cluster, counter].set_xlim([0, 70])
# 			axes[cluster,counter].set_xticks([0,70])
# 			counter += 1

# fig.tight_layout()
# plt.savefig("plots/smFISH_300_norm_med_traces.pdf")


"""
For a particular target, compute the mean and standard deviation for each cluster
"""
# target = "Ccl5"
# counts_dict = {"0": [], "1":[], "2":[]}

# for cell in good_cells:
# 	if cell.target == target:
# 		counts_dict[str(cell.clusterID)] += [cell.mRNA_count]

# for key in counts_dict.keys():
# 	counts_dict[key] = np.array(counts_dict[key])

# for key in ["0","1","2"]:
# 	print "Cluster " + key
# 	print "Mean: " + str(np.mean(counts_dict[key]))
# 	print "Std : " + str(np.std(counts_dict[key]))
# 	print "Std err: " + str(sem(counts_dict[key]))


"""
Visualize the mean mRNA counts with error bars over time
"""

times = ["0", "75", "150", "300"]
times_float = [0, 75, 150, 300]
targets = ["A20", "Cxcl2", "Ccl4", "Ccl5", "Nfkbia", "Nfkbiz", "Tnfsf9", "Atf3"]
# targets = ["A20", "Ccl5"]

plt.clf()
fig, axes = plt.subplots(len(targets)/2,2, figsize = (8*2,4*len(targets)))
counter = 0

master_dict = {}
for time in times:
	master_dict[time] = {}
	file_name_save = os.path.join(direc, "good_cells_" + time + "min.pkl")
	good_cells = pickle.load(open(file_name_save))
	for target in targets:
		counts_dict = {"0": [], "1":[], "2":[], "-1":[]}

		for cell in good_cells:
			if cell.target == target:
				value = cell.mRNA_count

				if cell.clusterID == None:
					if np.isnan(value) == False:
						counts_dict["0"] += [value]
						counts_dict["1"] += [value]
						counts_dict["2"] += [value]
				else:
					if np.isnan(value) == False:
						counts_dict[str(cell.clusterID)] += [value]


		for key in counts_dict.keys():
			counts_dict[key] = np.array(counts_dict[key])

		master_dict[time][target] = counts_dict

colors = {"0":"g", "1":"r", "2":"b"}
clusters = ["0", "1", "2"]


counter = 0
for target in targets:
	max_list_high = 0
	for cluster in clusters:
		list_of_times = []
		list_of_means = []
		list_of_errs = []
		for time in times:
				# Plot clusters 0, 1, 2, and -1
				points_to_plot = master_dict[time][target][cluster]
				if len(points_to_plot) > 0:
					mean = np.median(points_to_plot)
					std_err = 1.25*sem(points_to_plot)
					list_of_times += [np.float(time)]
					list_of_means += [mean]
					list_of_errs += [std_err]

		list_of_times = np.array(list_of_times)
		list_of_means = np.array(list_of_means)
		list_of_errs = np.array(list_of_errs)

		list_low = list_of_means - list_of_errs
		list_high = list_of_means + list_of_errs

		if np.amax(list_high) > max_list_high:
			max_list_high = np.amax(list_high)

		axes.flatten()[counter].plot(list_of_times, list_of_means, color = colors[cluster], linewidth = 1, label = "Cluster + " + cluster)
		axes.flatten()[counter].fill_between(list_of_times, list_low, list_high, color = colors[cluster], alpha = 0.1)
		axes.flatten()[counter].set_title(target, fontsize = 20)
		axes.flatten()[counter].set_xlim([0, 300])
		axes.flatten()[counter].set_xticks([0, 75, 150, 300])
		axes.flatten()[counter].set_ylim([0, max_list_high*1.05])
		axes.flatten()[counter].set_yticks([0, max_list_high*1.05])
		# axes.flatten()[counter].set_yticks([0, 100])


		axes.flatten()[counter].set_xlabel("Time (min)", fontsize = 20)
		axes.flatten()[counter].set_ylabel("Total cytoplasmic intensity (au)", fontsize = 20)
		axes.flatten()[counter].tick_params(axis = "both", which = "major", labelsize = 16)

		# axes.flatten()[counter].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
		axes.flatten()[counter].yaxis.major.formatter._useMathText = True
	counter += 1
fig.tight_layout()
# for key in ["0","1","2"]:
# 	print "Cluster " + key
# 	print "Mean: " + str(np.mean(master_dict["300"]["Ccl5"][key]))
# 	print "Std : " + str(np.std(master_dict["300"]["Ccl5"][key]))
# 	print "Std err: " + str(sem(master_dict["300"]["Ccl5"][key]))

plt.savefig("plots/smFISH_mRNA_counts.pdf")


"""
Plot correlation of different mRNA measures
"""
count = []
cyto = []
med = []
summed = []

for cell in good_cells:
	if np.isnan(cell.mRNA_count) == 0 and np.isnan(cell.mRNA_med_int) == 0:
		count += [cell.mRNA_count]
		cyto += [cell.cyto_fluo]
		med += [cell.mRNA_med_int]
		summed += [cell.mRNA_summed]

count = np.array(count)
cyto = np.array(cyto)
med = np.array(med)
summed = np.array(summed)

fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111)
slope, intercept, r_value, p_value, std_err = stats.linregress(count, cyto)
x = np.linspace(0,120,1000)
y = slope *x + intercept

ax.scatter(count, cyto, color = 'b', s = .5, alpha = 1)
ax.plot(x,y,'r')
ax.set_xlabel("mRNA counts")
ax.set_ylabel("Total cytoplasmic intensity")
fig.tight_layout()
ax.text(.3,.9, 'r = ' + str(r_value), ha='center', va='center', transform=ax.transAxes)
plt.savefig("plots/counts_vs_cyto.pdf")


plt.clf()
fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111)
slope, intercept, r_value, p_value, std_err = stats.linregress(count, summed)
x = np.linspace(0,120,1000)
y = slope *x + intercept

ax.scatter(count, summed, color = 'b', s = .5, alpha = 1)
ax.plot(x,y,'r')
ax.set_xlabel("mRNA counts")
ax.set_ylabel("Total mRNA intensity (segmented)")
ax.text(.3,.9, 'r = ' + str(r_value), ha='center', va='center', transform=ax.transAxes)
fig.tight_layout()
plt.savefig("plots/counts_vs_summed.pdf")

plt.clf()
fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111)
slope, intercept, r_value, p_value, std_err = stats.linregress(summed, cyto)
x = np.linspace(0,4e7,1000)
y = slope *x + intercept

ax.scatter(summed, cyto, color = 'b', s = .5, alpha = 1)
ax.plot(x,y,'r')
ax.set_xlabel("Total mRNA intensity (segmented)")
ax.set_ylabel("Total cytoplasmic intensity")
ax.text(.3,.9, 'r = ' + str(r_value), ha='center', va='center', transform=ax.transAxes)
fig.tight_layout()
plt.savefig("plots/summed_vs_cyto_fluo.pdf")

"""
Plot histograms of smFISH data points across clusters
"""

for target in targets:
	plt.clf()
	fig, ax = plt.subplots(3,1, figsize = (20,20), squeeze = False)
	for i in xrange(3):
		list_of_means = []
		list_of_errs = []
		cluster = clusters[i]


		for j in xrange(4):
			time = times[j]
			points_to_plot = master_dict[time][target][cluster]
			time_int = np.int(time)
			times_to_plot = np.ones(points_to_plot.shape)*time_int

			mean = np.mean(points_to_plot)
			err = sem(points_to_plot)
			list_of_means += [mean]
			list_of_errs += [err]

			ax[i,0].scatter(times_to_plot, points_to_plot, label = (time + " min"), color = colors[cluster])
			ax[i,0].set_title(target)
			ax[i,0].set_xlim([-5, 305])
			ax[i,0].set_ylim([-1e6, 3e7])
			ax[i,0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
		list_of_times = np.array([0, 75, 150, 300])
		ax[i,0].errorbar(list_of_times, list_of_means, yerr = list_of_errs, 
			fmt = "-o", color = colors[cluster], ecolor = colors[cluster], linewidth = 2, label = "Cluster + " + cluster)
	fig.tight_layout()
	plt.savefig("plots/smFISH/" + target + "_histograms.pdf")



"""
Plot heat map - normalized median
"""
# longest_time = 0
# number_of_cells = 0

# for cell in good_cells:
# 	number_of_cells += 1
# 	longest_time = np.amax([longest_time, cell.norm_med.shape[0]])

# heat_map = np.zeros((number_of_cells,longest_time))

# cell_counter = 0
# for cell in good_cells:
# 	if cell.clusterID == 0:
# 		number_of_cells += 1
# 		dynam = cell.norm_med
# 		heat_map[cell_counter,0:dynam.shape[0]] = dynam
# 		cell_counter += 1

# fig = plt.figure(figsize = (6,8))
# ax = fig.add_subplot(111)
# cleanAxis(ax)

# cax = ax.imshow(heat_map, cmap = plt.get_cmap('Reds'), interpolation = 'none')
# ax.set_xlabel('Time')
# ax.set_ylabel('Cells')
# ax.set_title('300 minute NFkB activity heatmap - smFISH cells - Norm med' + str(number_of_cells) + ' cells', y = 1.05)
# fig.colorbar(cax, ticks = [0, 1], orientation = 'vertical')

# plt.savefig("plots/smFISH_300_norm_med.pdf")

# """
# Load cells
# """
# direc = '/scratch/PI/mcovert/dvanva/sequencing/'
# all_cell_file = 'all_cells_all_detected_genes_qc_w_jackpot.pkl'
# all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

# """
# Select and cluster 300 min time point and plot some examples
# """

# time_point = 300

# print "Analyzing " + str(time_point) + " minute time point"
# all_cells = []
# longest_time = 0
# number_of_cells = 0

# for cell in all_cells_total:
# 	if cell.time_point == time_point and cell.condition == 'Stim':
# 		all_cells += [cell]
# 		number_of_cells += 1

# dynamics_load = np.load("/home/dvanva/SingleCellSequencing/300_dynamics_distance_matrix_kshape.npz")
# distance_matrix = dynamics_load['distance_matrix']
# Y_dynamics = sch.linkage(distance_matrix, method = 'ward')
# ind_dynamics = sch.fcluster(Y_dynamics,0.5*np.amax(Y_dynamics[:,2]),'distance')

# for j in xrange(number_of_cells):
# 	all_cells[j].clusterID = ind_dynamics[j]-1
# 	print all_cells[j].clusterID

# fig, axes = plt.subplots(3,10, figsize = (20,10))
# print axes.shape

# for cluster in [0,1,2]:
# 	counter = 0
# 	for cell in all_cells:
# 		if cell.clusterID == cluster and counter < 10:	
# 			print cluster, counter	
# 			axes[cluster, counter].plot(cell.NFkB_dynamics)
# 			axes[cluster, counter].set_xlim([0, 70])
# 			axes[cluster,counter].set_xticks([0,70])
# 			counter += 1

# fig.tight_layout()
# plt.savefig("plots/smFISH_300_c1_traces.pdf")
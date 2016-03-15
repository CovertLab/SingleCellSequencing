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
from seq_functions import smFISH_cell
import rpy2
from rpy2.robjects.packages import importr
import cPickle as pickle
rpy2.robjects.numpy2ri.activate()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42

"""
Initialize R instances
"""

R = rpy2.robjects.r
DTW = importr('dtw')
DTWCLUST = importr('dtwclust')

"""
Load excel files
"""

direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"

file_name = os.path.join(direc, "12072015", "12072015_ExperMetadata.xlsx")
data_0 = pd.read_excel(file_name, sheetname = 0)

file_name = os.path.join(direc, "12162015", "12162015_ExperMetadata.xlsx")
data_1 = pd.read_excel(file_name, sheetname = 0)

"""
Load cluster averages
"""
direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name = os.path.join(direc,"300_cluster_avg_kshape_c1.npz")
file_load = np.load(file_name)

"""
Load MAT files
"""

# times = ["0", "75", "150", "300"]
times = ["300"]
dates = ["12072015", "12162015"]
for time in times:
	print time

	if time == "0":
		cluster_dynamics_avg = None

	else:
		cluster_dynamics_avg = file_load["cluster_dynamics_avg"]
		file_name = os.path.join(direc, "12072015", "nucData_withFISH_3.mat")
		dynamics_file = sio.loadmat(file_name)
		temp = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,7:]
		longest_time = temp.shape[1]
		cluster_dynamics_avg = cluster_dynamics_avg[:,0:longest_time]
	
	list_of_cells = []

	for date in dates:
		file_name = os.path.join(direc, date, "nucData_withFISH_3.mat")
		dynamics_file = sio.loadmat(file_name)

		conditions = dynamics_file["dataToAnalyzeNuc"+ time + "Ratio"][:,0]
		positions = dynamics_file["dataToAnalyzeNuc"+ time + "Ratio"][:,1]
		ids = dynamics_file["dataToAnalyzeNuc"+ time + "Ratio"][:,2]

		counts = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,3]
		medint = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,4]
		summed = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,5]
		cyto_fluo = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,6]

		ratios = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,7:]
		normalized = dynamics_file["normalizedNuc" + time + "Data"][:,7:]

		print conditions.shape
		number_of_cells = conditions.shape[0]

		for j in xrange(number_of_cells):
			if date == "12072015":
				cell = smFISH_cell(cell_id = ids[j], NC_ratio = ratios[j,:], norm_med = normalized[j,:], condition = conditions[j], counts = counts[j], 
					medint = medint[j], cyto_fluo = cyto_fluo[j], summed = summed[j], position = positions[j], smFISH_dataframe = data_0, cluster_dynamics_avg = cluster_dynamics_avg)

			if date == "12162015":
				cell = smFISH_cell(cell_id = ids[j], NC_ratio = ratios[j,:], norm_med = normalized[j,:], condition = conditions[j], counts = counts[j], 
					medint = medint[j], cyto_fluo = cyto_fluo[j], summed = summed[j], position = positions[j], smFISH_dataframe = data_1, cluster_dynamics_avg = cluster_dynamics_avg)

			list_of_cells += [cell]

	good_cells = []
	for cell in list_of_cells:
		if cell.good_cell == 1:
			good_cells += [cell]

	# if time == "300":
	# 	dynamics_load = np.load("/scratch/PI/mcovert/dvanva/sequencing/smFISH/300_dynamics_distance_matrix_kshape.npz")
	# 	distance_matrix = dynamics_load['distance_matrix']
	# 	Y = sch.linkage(distance_matrix, method = 'ward')
	# 	ind_dynamics = sch.fcluster(Y,0.5*np.amax(Y[:,2]),'distance') - 1

	# 	for j in xrange(len(good_cells)):
	# 		good_cells[j].clusterID = ind_dynamics[j]


	file_name_save = os.path.join(direc, "good_cells_" + time + "min.pkl")
	pickle.dump(good_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

	print len(list_of_cells), len(good_cells)

	"""
	Plot heat map using the c1 clustering
	"""

	"""
	Fill up the heat map matrix
	"""

	longest_time = 0
	number_of_cells = 0
	for cell in good_cells:
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.norm_med.shape[0]])

	dynamics_matrix = np.zeros((number_of_cells,longest_time))

	if time in ["75", "150", "300"]:
		cluster_len_1 = 0
		cluster_len_2 = 0
		cluster_len_3 = 0

		for cell in good_cells:
			if cell.clusterID == 0:
				cluster_len_1 += 1
			if cell.clusterID == 1:
				cluster_len_2 += 1
			if cell.clusterID == 2:
				cluster_len_3 += 1

		print cluster_len_1, cluster_len_2, cluster_len_3
		frac_1 = np.float(cluster_len_1) / np.float(cluster_len_1 + cluster_len_2 + cluster_len_3)
		frac_2 = np.float(cluster_len_2) / np.float(cluster_len_1 + cluster_len_2 + cluster_len_3)
		frac_3 = np.float(cluster_len_3) / np.float(cluster_len_1 + cluster_len_2 + cluster_len_3)

		print frac_1, frac_2, frac_3

		cell_counter = 0

		fig = plt.figure(figsize = (8,4))
		ax1 = fig.add_subplot(131)
		ax2 = fig.add_subplot(132)
		ax3 = fig.add_subplot(133)

		counter_max = 1
		times = np.arange(0,longest_time)*5

		colors = ['g', 'r', 'b','purple','yellow']

		counter = 0
		for cell in good_cells:
			if cell.clusterID == 0: 
				if counter == 5:
					ax3.plot(times, cell.norm_med, color = colors[0])
				counter += 1
		
		counter = 0
		for cell in good_cells:
			if cell.clusterID == 1:
				if counter == 6:
					ax2.plot(times, cell.norm_med, color = colors[1])
				counter += 1
		
		counter = 0
		for cell in good_cells:
			if cell.clusterID == 2:
				if counter == 7:
					ax1.plot(times, cell.norm_med, color = colors[2])
				counter += 1

		ax1.set_ylabel('Nuclear localization (au)', fontsize = 16)
		ax1.set_xlim([0, times[-1]])
		ax2.set_xlim([0,times[-1]])
		ax3.set_xlim([0,times[-1]])

		ax1.set_ylim([0,1])
		ax2.set_ylim([0,1])
		ax3.set_ylim([0,1])

		ax1.set_xticks([0, 300])
		ax2.set_xticks([0, 300])
		ax3.set_xticks([0, 300])

		ax1.set_yticks([0,1])
		ax2.set_yticks([0,1])
		ax3.set_yticks([0,1])

		plt.tight_layout()

		# ax_heatmap_1 = fig.add_axes([0.3, 0.1 +0.005, 0.6, 0.8*frac_1 - 0.005], frame_on = True)
		# ax_heatmap_2 = fig.add_axes([0.3, 0.1 +0.005+ 0.8*frac_1, 0.6 , 0.8*frac_2 - 0.005], frame_on = True)
		# ax_heatmap_3 = fig.add_axes([0.3, 0.1 +0.005+ 0.8*frac_1 + 0.8*frac_2, 0.6, 0.8*frac_3 - 0.005], frame_on = True)

		# for cell in good_cells:
		# 	if cell.clusterID == 0:
		# 		dynam = cell.norm_med
		# 		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		# 		cell_counter += 1

		# for cell in good_cells:
		# 	if cell.clusterID == 1:
		# 		dynam = cell.norm_med
		# 		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		# 		cell_counter += 1

		# for cell in good_cells:
		# 	if cell.clusterID == 2:
		# 		dynam = cell.norm_med
		# 		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		# 		cell_counter += 1


		# # ax_heatmap = fig.add_axes([0.3, 0.1, 0.6, 0.8])
		# im1 = ax_heatmap_1.matshow(dynamics_matrix[0:cluster_len_1,:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')
		# im2 = ax_heatmap_2.matshow(dynamics_matrix[cluster_len_1:cluster_len_1+cluster_len_2,:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')
		# im3 = ax_heatmap_3.matshow(dynamics_matrix[cluster_len_1+cluster_len_2:,:], aspect = 'auto', origin = 'lower', cmap = plt.get_cmap('Reds'), interpolation = 'none')


		# # fig.colorbar(im1, ticks = [0, 1], orientation = 'vertical')

		# # ax_heatmap.xaxis.set_ticks(np.arange(63))
		# ax_heatmap_1.xaxis.set_ticklabels([])
		# ax_heatmap_1.xaxis.set_ticks([])
		# ax_heatmap_1.yaxis.set_ticklabels([])
		# ax_heatmap_1.yaxis.set_ticks([])
		# ax_heatmap_2.yaxis.set_ticklabels([])
		# ax_heatmap_2.yaxis.set_ticks([])
		# ax_heatmap_2.xaxis.set_ticklabels([])
		# ax_heatmap_2.xaxis.set_ticks([])
		# ax_heatmap_3.yaxis.set_ticklabels([])
		# ax_heatmap_3.yaxis.set_ticks([])
		# ax_heatmap_3.xaxis.set_ticklabels([])
		# ax_heatmap_3.xaxis.set_ticks([])
		# ax_heatmap_3.set_title('smFISH - ' + str(time) + " minutes - " + str(number_of_cells) + ' cells', y = 1.05)
		# ax_heatmap_1.set_xlabel('Time (minutes)')

		# plt.savefig("plots/" + time + "min_smFISH_dynamics_c1_clustering.pdf")
		plt.savefig("plots/" + time + "min_smFISH_dynamics_c1_clustering_samples.pdf")





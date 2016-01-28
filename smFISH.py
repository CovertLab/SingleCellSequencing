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

"""
Initialize R instances
"""

R = rpy2.robjects.r
DTW = importr('dtw')
DTWCLUST = importr('dtwclust')

"""
Load cluster averages
"""
direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name = os.path.join(direc,"300_cluster_avg_kshape_smFISH.npz")
file_load = np.load(file_name)

"""
Load excel files
"""

direc = "/scratch/PI/mcovert/dvanva/sequencing/smFISH"
file_name = os.path.join(direc,"smFISH_counts.xlsx")

data_0 = pd.read_excel(file_name, sheetname = 0)
data_1 = pd.read_excel(file_name, sheetname = 1)

"""
Load MAT files
"""

times = ["0", "75", "150", "300"]
dates = ["12072015", "12162015"]
for time in times:
	print time

	if time == "0":
		cluster_dynamics_avg = None

	else:
		cluster_dynamics_avg = file_load["cluster_dynamics_avg"]
		file_name = os.path.join(direc, "12072015", "nucData.mat")
		dynamics_file = sio.loadmat(file_name)
		temp = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,3:]
		longest_time = temp.shape[1]
		cluster_dynamics_avg = cluster_dynamics_avg[:,0:longest_time]
	
	list_of_cells = []

	for date in dates:
		file_name = os.path.join(direc, date, "nucData.mat")
		dynamics_file = sio.loadmat(file_name)

		conditions = dynamics_file["dataToAnalyzeNuc"+ time + "Ratio"][:,0]
		positions = dynamics_file["dataToAnalyzeNuc"+ time + "Ratio"][:,1]
		ids = dynamics_file["dataToAnalyzeNuc"+ time + "Ratio"][:,2]

		ratios = dynamics_file["dataToAnalyzeNuc" + time + "Ratio"][:,3:]
		normalized = dynamics_file["normalizedNuc" + time + "Data"][:,3:]

		number_of_cells = conditions.shape[0]

		for j in xrange(number_of_cells):
			if date == "12072015":
				cell = smFISH_cell(cell_id = ids[j], NC_ratio = ratios[j,:], norm_med = normalized[j,:], condition = conditions[j], 
					position = positions[j], smFISH_dataframe = data_0, cluster_dynamics_avg = cluster_dynamics_avg)

			if date == "12162015":
				cell = smFISH_cell(cell_id = ids[j], NC_ratio = ratios[j,:], norm_med = normalized[j,:], condition = conditions[j], 
					position = positions[j], smFISH_dataframe = data_1, cluster_dynamics_avg = cluster_dynamics_avg)

			list_of_cells += [cell]

	good_cells = []
	for cell in list_of_cells:
		if cell.good_cell == 1:
			good_cells += [cell]

	file_name_save = os.path.join(direc, "good_cells_" + time + "min.pkl")
	pickle.dump(good_cells, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)

	print len(good_cells)


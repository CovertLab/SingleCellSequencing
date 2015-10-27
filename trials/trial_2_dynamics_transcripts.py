"""
Analysis!

Plot the dynamics for the first 3 cells in the all_cell object
Also plot the TPMs for the most abundant transcripts (10) in each cell's genome
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
import matplotlib
import cPickle as pickle
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyensembl

"""
Load all the cells
"""

matplotlib.style.use('ggplot')

# mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')


direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

cell_counter = 0 
num_cells = 8

fig, axes = plt.subplots(2,num_cells, figsize = (28,10))

for cell in all_cells:
	if cell_counter > num_cells - 1:
		break

	if cell.time_point == 150:
		axes[0,cell_counter].plot(cell.NFkB_dynamics, label = 'Cell ' + str(cell_counter + 1))
		axes[0,cell_counter].set_xlabel('Time (time points)')
		axes[0,cell_counter].set_ylabel('Normalized fluorescence (au)')
		axes[0,cell_counter].set_title('NFkB dynamics for cell ' + str(cell_counter + 1), fontsize = 12)

		sorted_transcripts = cell.transcriptome.sort(columns = 'tpm', ascending = False)
		top_ten_tpm = sorted_transcripts[0:10].tpm
		top_ten_indices = sorted_transcripts[0:10].index.tolist()

		# for i in xrange(len(top_ten_indices)):
		# 	top_ten_indices[i] = mouse_genome.gene_name_of_transcript_id(top_ten_indices[i])

		ind = np.arange(10)
		width = 0.35
		axes[1,cell_counter].bar(ind, top_ten_tpm, width)
		axes[1,cell_counter].set_ylabel('TPM (transcripts per million)')
		axes[1,cell_counter].set_title('TPMs for cell ' + str(cell_counter + 1), fontsize = 12)

		axes[1,cell_counter].set_xticks(ind + width/2)
		xtickNames = axes[1,cell_counter].set_xticklabels(top_ten_indices)
		plt.setp(xtickNames, rotation = 45, fontsize = 8, ha = 'right')

		cell_counter += 1
fig.tight_layout()

plt.savefig("trial_2_dynamics_transcripts_with_gene_names.pdf")




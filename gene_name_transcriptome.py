"""
Construct gene name table
The transcriptome table is currently a list of estimated counts for each transcript.
However, some genes have multiple isoforms that are represented as separate transcripts.
This script constructs a new transcriptome table that can be addressed by gene name by
summing up the counts of each gene over all of its isoforms.

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
import seq_functions

matplotlib.style.use('ggplot')

"""
Load cells
"""

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc.pkl'
all_cells = pickle.load(open(os.path.join(direc,all_cell_file)))

"""
Load mouse genome
"""

mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')
gene_names = pyensembl.gene_names

transcript_dict = {}

for gene in gene_names:
	transcript_dict[gene] = transcript_ids_of_gene_name(gene)

for cell in all_cells:
	counts_dict = {}
	tpm_dict = {}
	for gene in gene_names:
		transcript_names = transcript_dict[gene]
		counts = 0
		tpm = 0
		for transcript_name in transcript_names:
			counts += cell.transcripts.loc[transcript_name]['est_counts']
			counts += cell.transcripts.loc[transcript_name]['tpm']


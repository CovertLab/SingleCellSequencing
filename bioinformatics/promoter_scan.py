"""
Analysis!

Pull out the promoters from a list of genes
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
from pyfaidx import Fasta as fasta
import Bio
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Alphabet import IUPAC
import matplotlib.ticker as mtick


"""
Load mouse genome
"""

mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')

# list_of_genes = ["Nfkbia", "Tnfaip3", "Nfkbiz", "Nfkbie", "Tnfaip2", "Nlrp3", "Ccl3", "Ccl5", "Saa3", "Il1f9", "Il6", "Csf3", "Gp49a", "Cxcl2", "Cxcl10", "Ccl20", "Tnfsf9", "Il1a", "Tnf", "Ccl4", "Il1b", "Cxcl3", "Il1f6", "Lif"]

list_of_lists = []
list_of_genes = ["Nfkbia", "Tnfaip3", "Nfkbiz", "Nfkbie", "Tnfaip2", "Nlrp3"]
list_of_lists += [list_of_genes]
list_of_genes = ["Ccl3", "Ccl5"]
list_of_lists += [list_of_genes]

list_of_genes = ["Saa3", "Il1f9", "Il6", "Csf3", "Gp49a"]
list_of_lists += [list_of_genes]

list_of_genes = ["Cxcl2", "Cxcl10", "Ccl20", "Tnfsf9", "Il1a"]
list_of_lists += [list_of_genes]

list_of_genes = ["Tnf", "Ccl4", "Il1b", "Cxcl3", "Il1f6", "Lif"]
list_of_lists += [list_of_genes]

names = ["feedback", "cluster2", "osc", "75min", "150min"]
list_of_tfs = ["RELA", "REL", "NFKB2", "Ddit3::Cebpa", "FOS::JUN", "FOS", "JUN", "STAT1", "STAT1::STAT2", "STAT3", "Stat4", "Stat5a::Stat5b", "Stat6", "Foxo1", "FOXO3", "Bcl6", "Atf3", "CEBPA", "CEBPB", "CEBPD", "NFAT5", "NFATC1", "CREB1"]
transcript_dict = {}
gene_object_dict = {}
promoter_dict = {}

for list_of_genes, name in zip(list_of_lists, names):
	print name
	print list_of_genes
	for gene in list_of_genes:
		transcript_dict[gene] = mouse_genome.transcript_ids_of_gene_name(gene)

	for gene in list_of_genes:
		gene_object_dict[gene] = mouse_genome.genes_by_name(gene)[0]

	"""
	Look up promoter sequences
	"""

	genome_loc = "/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/mouse_genome_w_spikeins.fa"
	chromosomes = fasta(genome_loc)

	upstream = 1000
	downstream = 500
	for gene_name in list_of_genes:
		gene = gene_object_dict[gene_name]
		if gene.strand == "+":
			start_pos = gene.start-1
			promoter = Seq(chromosomes[str(gene.contig)][start_pos - upstream:start_pos + downstream].seq, IUPAC.unambiguous_dna)
		elif gene.strand == "-":
			start_pos = gene.end
			promoter = Seq(chromosomes[str(gene.contig)][start_pos - downstream:start_pos + upstream].reverse.complement.seq, IUPAC.unambiguous_dna)
		promoter_dict[gene_name] = promoter

	"""
	Load JASPAR motifs
	"""

	motif_dict = {}
	ps = {"A": 0.6, "C": 0.4, "T": 0.6, "G": 0.4}
	fh = open("pfm_vertebrates.txt")
	for m in motifs.parse(fh, "jaspar"):
		pwm = m.counts.normalize(pseudocounts = ps)
		pssm = pwm.log_odds()
		motif_dict[m.name[1:]] = pssm 

	list_of_motifs = motif_dict.keys()

	# print motif_dict["NFKB1"].calculate(p_seq)

	plt.clf()
	fig, axes = plt.subplots(len(list_of_genes),len(list_of_tfs), figsize = (5*len(list_of_tfs),5*len(list_of_genes)))
	counter_2 = 0
	for transcription_factor in list_of_tfs:
		counter_1 = 0

		print transcription_factor
	# transcription_factor = "RELA"
		hits_dict = {}
		hits_search_dict = {}
		threshold = 8

	# print motif_dict[transcription_factor].consensus

		for gene in list_of_genes:
			hits_dict[gene] = []

		for gene in list_of_genes:
			motif_length = len(motif_dict[transcription_factor].consensus)
			prediction = motif_dict[transcription_factor].calculate(promoter_dict[gene])
			hits_search_dict[gene] = motif_dict[transcription_factor].search(promoter_dict[gene], threshold = threshold)

			for pos, score in hits_search_dict[gene]:
				hits_dict[gene]	+= [promoter_dict[gene][pos:pos + motif_length]]
			position = np.arange(-upstream,downstream)[0:-len(motif_dict[transcription_factor].consensus)+1]

			final_pred = 2** prediction / 2**motif_dict[transcription_factor].max
			max_val = np.amax(final_pred[~np.isnan(final_pred)])
			if max_val > .01:
				axes[counter_1, counter_2].plot(position, final_pred , "g", linewidth = 1)
			axes[counter_1, counter_2].set_title(gene + " , " + transcription_factor)
			axes[counter_1,counter_2].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
			counter_1 +=1
		counter_2 +=1
		# print hits_dict
	fig.tight_layout()

	file_name = "prom_scan_" + name + ".pdf"
	plt.savefig("plots/" + file_name)



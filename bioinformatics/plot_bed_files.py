"""
Analysis!

Plot bed graph tracks
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
import track

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42

"""
Convert tracks to sql (only has to be done once)
"""

# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM487448_p300.UT.peaks.bed','/scratch/PI/mcovert/dvanva/bedfiles/sql/p300.UT.peaks.sql')
# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM487449_p300.LPS_2h.peaks.bed', '/scratch/PI/mcovert/dvanva/bedfiles/sql/p300.LPS_2h.peaks.sql')
# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM487450_PU.1.UT.peaks.bed', '/scratch/PI/mcovert/dvanva/bedfiles/sql/PU.1.UT.peaks.sql' )
# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM487451_PU.1.LPS_2h.peaks.bed', '/scratch/PI/mcovert/dvanva/bedfiles/sql/PU.1.LPS_2h.peaks.sql')
# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM487452_H3K4me1.UT.peaks.bed', '/scratch/PI/mcovert/dvanva/bedfiles/sql/H3K4me.UT.peaks.sql')

# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM1913601_0m_rep3.bedGraph', '/scratch/PI/mcovert/dvanva/bedfiles/sql/ATAC_0m_rep3.sql')
# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM1913605_120m_rep3.bedGraph', '/scratch/PI/mcovert/dvanva/bedfiles/sql/ATAC_120m_rep3.sql')

# track.convert('/scratch/PI/mcovert/dvanva/bedfiles/GSM1645121_RelA-2_120.bedgraph', '/scratch/PI/mcovert/dvanva/bedfiles/sql/RelA_120m_rep2.sql')

"""
Plot tracks around Ccl3 and Ccl4
"""

plt.clf()
fig, axes = plt.subplots(10,1, figsize = (6,1*10))
start_pos = 83460000
end_pos = 83480000

# p300 UT peaks
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/p300.UT.peaks.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[0].bar(peak[0], 1, width = peak[1]-peak[0], color = 'k')
		axes[0].set_xlim([start_pos, end_pos])
		axes[0].set_title('p300 Untreated')
		axes[0].set_xticks([])
		axes[0].set_yticks([])

# p300 LPS 2h peaks
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/p300.LPS_2h.peaks.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[1].bar(peak[0], 1, width = peak[1]-peak[0], color = 'k')
		axes[1].set_xlim([start_pos, end_pos])
		axes[1].set_title('p300 LPS 2h')
		axes[1].set_xticks([])
		axes[1].set_yticks([])

# PU.1 UT peaks
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/PU.1.UT.peaks.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[2].bar(peak[0], 1, width = peak[1]-peak[0], color = 'k')
		axes[2].set_xlim([start_pos, end_pos])
		axes[2].set_title('PU.1 Untreated')
		axes[2].set_xticks([])
		axes[2].set_yticks([])


# PU.1 LPS 2h peaks
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/PU.1.LPS_2h.peaks.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[3].bar(peak[0], 1, width = peak[1]-peak[0], color = 'k')
		axes[3].set_xlim([start_pos, end_pos])
		axes[3].set_title('PU.1 LPS 2h')
		axes[3].set_xticks([])
		axes[3].set_yticks([])


# H3K4me1 peaks
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/H3K4me.UT.peaks.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[4].bar(peak[0], 1, width = peak[1]-peak[0], color = 'k')
		axes[4].set_xlim([start_pos, end_pos])
		axes[4].set_title('H3K4me1 Untreated')
		axes[4].set_xticks([])
		axes[4].set_yticks([])


# ATAC seq 0 min
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/ATAC_0m_rep3.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[6].bar(peak[0], peak[2], width = peak[1]-peak[0], color = 'b', edgecolor = 'none')
		axes[6].set_xlim([start_pos, end_pos])
		axes[6].set_title('ATAC-Seq Untreated')
		axes[6].set_xticks([])
		axes[6].set_yticks([0, 12])


# ATAC seq 120 min
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/ATAC_120m_rep3.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[7].bar(peak[0], peak[2], width = peak[1]-peak[0], color = 'b', edgecolor = 'none')
		axes[7].set_xlim([start_pos, end_pos])
		axes[7].set_title('ATAC-Seq Lipid A 2h')
		axes[7].set_xticks([])
		axes[7].set_yticks([0, 50])


# RelA  120 min
with track.load('/scratch/PI/mcovert/dvanva/bedfiles/sql/RelA_120m_rep2.sql') as t:
	data = t.read({'chr':'chr11', 'start': 83460000, 'end': 83480000})
	for peak in data:
		axes[8].bar(peak[0], peak[2], width = peak[1]-peak[0], color = 'b', edgecolor = 'none')
		axes[8].set_xlim([start_pos, end_pos])
		axes[8].set_title('RelA Lipid A 2h')
		axes[8].set_yticks([0, 40])
		axes[8].set_xticks([])

"""
Plot location of Ccl3 and Ccl4
"""
mouse_genome = pyensembl.Genome(reference_name = 'NCBIM37', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/cdna/Mus_musculus.NCBIM37.67.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/pep/Mus_musculus.NCBIM37.67.pep.all.fa.gz')
list_of_genes = ['Ccl3', 'Ccl4']
transcript_dict = {}
gene_object_dict = {}

for gene in list_of_genes:
	transcript_dict[gene] = mouse_genome.transcript_ids_of_gene_name(gene)

for gene in list_of_genes:
	gene_object_dict[gene] = mouse_genome.genes_by_name(gene)[0]

# Gene locations
for gene_name in list_of_genes:
	gene = gene_object_dict[gene_name]
	gene_start = gene.start
	width = gene.end - gene.start

	axes[5].bar(gene_start, 1, width = width, color = 'g', edgecolor = 'none')
	axes[5].set_xlim([start_pos, end_pos])
	axes[5].set_xticks([])
	axes[5].set_yticks([])
	axes[5].set_title('Coding regions')

"""
Look up promoter sequences
"""

genome_loc = "/home/dvanva/Mus_musculus.NCBIM37.67.dna.chromosome.11.fa"
chromosomes = fasta(genome_loc)
ccl3 = gene_object_dict['Ccl3']
intergenic = Seq(chromosomes[str(ccl3.contig)][start_pos:end_pos].seq, IUPAC.unambiguous_dna)
# print intergenic

motif_dict = {}
ps = {"A": 0.6, "C": 0.4, "T": 0.6, "G": 0.4}
fh = open("pfm_vertebrates.txt")
for m in motifs.parse(fh, "jaspar"):
	pwm = m.counts.normalize(pseudocounts = ps)
	pssm = pwm.log_odds()
	motif_dict[m.name[1:]] = pssm 

list_of_motifs = motif_dict.keys()

prediction = 2**motif_dict["RELA"].calculate(intergenic)
minimum = 2** motif_dict["RELA"].min
maximum = 2**motif_dict["RELA"].max
final_pred = (prediction - minimum)/(maximum-minimum)

# print prediction.shape
# print np.arange(start_pos, end_pos)[4:-5].shape

axes[9].plot(np.arange(start_pos, end_pos)[4:-5], final_pred, color = 'c')

prediction = 2**motif_dict["RELA"].calculate(intergenic.reverse_complement())
minimum = 2** motif_dict["RELA"].min
maximum = 2**motif_dict["RELA"].max
final_pred = (prediction - minimum)/(maximum-minimum)

axes[9].plot(np.arange(start_pos, end_pos)[4:-5], -final_pred, color = 'c')
axes[9].set_xlim([start_pos, end_pos])
axes[9].set_xticks([start_pos, end_pos])
axes[9].set_ylim([-0.04, 0.04])
axes[9].set_yticks([-0.04, 0.04])

axes[9].set_title('RelA PSSM Score')


for position, score in motif_dict["RELA"].search(intergenic, threshold = 8.0):
	print position, score

print intergenic[12492:12492+10]
fig.tight_layout()
file_name = "bed_plots.pdf"
plt.savefig("plots/" + file_name)


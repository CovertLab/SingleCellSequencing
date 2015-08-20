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

import seq_functions

"""
Define directories
"""

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
direc_mat = direc + 'imaging/'
library = 'library1/'

"""
Make directories
"""
# make_directories(direc + library)

"""
Unzip directory
"""
# unzip_rawdata(direc + library)

"""
Trim reads
"""

# trim_direc(direc + library)

"""
Make kallisto index file
"""
# direc_ref = 'reference_sequences/'
# mouse_fasta = 'Mus_musculus.GRCm38.rel79.cdna.all.fa.gz'
# input_filename = direc + direc_ref + mouse_fasta
# output_filename = direc + direc_ref + 'mouse'
# run_kallisto_test()

"""
Test sorting sam files
"""
# direc_name =  '/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/Mus_musculus.GRCm38.81.gtf'
# seq_functions.load_ensembl_gene_ids(direc_name)

file_name = '/scratch/PI/mcovert/dvanva/seq_test_run/library1/'

transcripts, spikes = seq_functions.count_hdf5_files(file_name, spikeids = ['Spike1', 'Spike4', 'Spike7'])
print spikes
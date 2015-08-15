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

from seq_functions import run_cmd, write_file, parse_filename, make_directories, unzip_file, unzip_rawdata, load_matfiles
from seq_functions import trim_direc, generate_genome_index, run_STAR, load_sequence_counts, run_kallisto_test

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
direc_ref = 'reference_sequences/'
mouse_fasta = 'Mus_musculus.GRCm38.rel79.cdna.all.fa.gz'
input_filename = direc + direc_ref + mouse_fasta
output_filename = direc + direc_ref + 'mouse'
run_kallisto_test()





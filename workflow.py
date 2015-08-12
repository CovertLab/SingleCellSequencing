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

from seq_functions import run_cmd, write_file, parse_filename, make_directories, unzip_file, unzip_rawdata
from seq_functions import load_bbmap, load_STAR, trim_direc, generate_genome_index, run_STAR, load_sequence_counts

"""
Define directories
"""

direc = '/scratch/users/dvanva/Sequencing/'
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
Create genome index
"""
generate_genome_index()


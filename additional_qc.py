"""
Additional qc - remove low expression remove_low_expression_genes
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
import cPickle as pickle
import seq_functions

"""
Load cells
"""

direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc.pkl'
all_cells_qc = pickle.load(open(os.path.join(direc,all_cell_file)))

all_cells_qc = seq_functions.remove_low_expression_genes(all_cells_qc)

file_name_save = os.path.join(direc, 'all_cells_qc_complete.pkl')
pickle.dump(all_cells_qc, open(file_name_save, 'wb'), protocol = pickle.HIGHEST_PROTOCOL)
"""
Analysis!

Load the listsdifferentially expressed genes in pandas
"""

"""
Import python packages
"""

import os
import pandas as pd
import numpy as np
import matplotlib
import cPickle as pickle
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

"""
Load scde results into data frames
"""

direc = "/scratch/PI/mcovert/dvanva/sequencing/"

file_1 = os.path.join(direc,"300min_scde_results_1v2.txt")
file_2 = os.path.join(direc,"300min_scde_results_1v3.txt")

gene_list_1v2 = pd.read_table(file_1, sep = "\t")
gene_list_1v3 = pd.read_table(file_2, sep = "\t")

"""
Sort by abs(Z) and then by abs(mle)
"""

new_index = gene_list_1v2.sort_values(["mle","Z"], ascending = [False, False], inplace = False).index
gene_list_1v2 = gene_list_1v2.reindex(index = new_index)

new_index = gene_list_1v3.sort_values(["mle","Z"], ascending = [False, False], inplace = False).index
gene_list_1v3 = gene_list_1v3.reindex(index = new_index)

trunc_1v2 = gene_list_1v2.iloc[0:200]
trunc_1v3 = gene_list_1v3.iloc[0:200]

print trunc_1v2.index
print trunc_1v3.index

set_1v2 = set(trunc_1v2.index)
set_1v3 = set(trunc_1v3.index)

final_set = set_1v2
final_set &= set_1v3

print list(final_set)


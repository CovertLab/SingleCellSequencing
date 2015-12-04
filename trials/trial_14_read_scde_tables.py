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

file_75 = os.path.join(direc,"75min_scde_results.txt")
file_150 = os.path.join(direc,"150min_scde_results.txt")
file_300 = os.path.join(direc,"300min_scde_results.txt")

gene_list_75 = pd.read_table(file_75, sep = "\t")
gene_list_150 = pd.read_table(file_150, sep = "\t")
gene_list_300 = pd.read_table(file_300, sep = "\t")

"""
Sort by abs(Z) and then by abs(mle)
"""

new_index = gene_list_75.abs().sort_values(["mle","Z"], ascending = [False, False], inplace = False).index
gene_list_75 = gene_list_75.reindex(index = new_index)

new_index = gene_list_150.abs().sort_values(["mle","Z"], ascending = [False, False], inplace = False).index
gene_list_150 = gene_list_150.reindex(index = new_index)

new_index = gene_list_300.abs().sort_values(["mle","Z"], ascending = [False, False], inplace = False).index
gene_list_300 = gene_list_300.reindex(index = new_index)

trunc_75 = gene_list_75.iloc[0:100]
trunc_150 = gene_list_150.iloc[0:100]
trunc_300 = gene_list_300.iloc[0:100]

print trunc_75.index
print trunc_150.index
print trunc_300.index

set_75 = set(trunc_75.index)
set_150 = set(trunc_150.index)
set_300 = set(trunc_300.index)

final_set = set_75
final_set &= set_150
final_set &= set_300

print final_set 

inflammatory_genes = ["Cxcl3", "Cxcl2", "Lif", "Ccl4", "Csf3", "Il1f9", "Ccl3", "Ccl5", "Tnf", "Il1a", "Il1b", "Tnfsf9", "Ccl20", "Il1f6", "Il27", "Il6"]
regulatory_genes = ["Nlrp3", "Nfkbiz", "Tnfaip2", "Nfkbia", "Tnfaip3", "Nfatc1"]
metabolic_genes = ["Hmox", "Prdx1", "Hdc", "Ptgs2", "Irg1"]
other_genes = ["Plaur", "Sqstm1", "Clec4e", "Sdc4", "Procr", "Slpi", "Plk2", "Saa3", "Slc7a11", "Cish", "Gp49a", "Hcar2", "Gpr84", "Malt1"]
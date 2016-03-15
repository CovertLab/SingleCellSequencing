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

new_index = gene_list_75.sort_values(["Z","mle"], ascending = [True, True], inplace = False).index
gene_list_75 = gene_list_75.reindex(index = new_index)

new_index = gene_list_150.sort_values(["Z","mle"], ascending = [True, True], inplace = False).index
gene_list_150 = gene_list_150.reindex(index = new_index)

new_index = gene_list_300.sort_values(["Z","mle"], ascending = [True, True], inplace = False).index
gene_list_300 = gene_list_300.reindex(index = new_index)

trunc_75 = gene_list_75.iloc[0:75]
trunc_150 = gene_list_150.iloc[0:75]
trunc_300 = gene_list_300.iloc[0:75]

set_75 = set(trunc_75.index)
set_150 = set(trunc_150.index)
set_300 = set(trunc_300.index)


# print set_75
# print set_150
# print set_300

print sorted(list(set_75 | set_150 | set_300))
print len(set_75 | set_150 | set_300)

final_set = set_75
final_set &= set_150
final_set &= set_300

# print final_set 

inflammatory_genes = ["Cxcl3", "Cxcl2", "Lif", "Ccl4", "Csf3", "Il1f9", "Ccl3", "Ccl5", "Tnf", "Il1a", "Il1b", "Tnfsf9", "Ccl20", "Il1f6", "Il27", "Il6"]
regulatory_genes = ["Nlrp3", "Nfkbiz", "Tnfaip2", "Nfkbia", "Tnfaip3", "Nfatc1"]
metabolic_genes = ["Hmox", "Prdx1", "Hdc", "Ptgs2", "Irg1"]
other_genes = ["Plaur", "Sqstm1", "Clec4e", "Sdc4", "Procr", "Slpi", "Plk2", "Saa3", "Slc7a11", "Cish", "Gp49a", "Hcar2", "Gpr84", "Malt1"]

all_genes = ['Wrb', 'Cxcl3', 'Tnfaip3', 'Tnfaip2', 'Gm6377', 'Blvrb', 'Bcl6b', 'Gpr21', 'Rsl1', 'Nfatc1', 'Nfkbiz', 'Phf11d', 'Tfec', 'Phf11a', 'Fas', 'Gpx1', 'Ccnd1', 'Flrt3', 'Cd247', 'Tslp', 'Atp6v0d2', 'Il1f6', 'Rnf213', 'Marcksl1', 'Abcg1', 'Gm8818', 'Nlrp3', 'Pde6b', 'Rrm2', 'Mbnl1', 'Ptpn14', 'Odc1', 'Fv1', 'Ptgs2', 'Ddit3', 'Slc25a2', 'Cfb', 'Mmp3', 'Ptprg', 'Cxcl10', 'Ly86', 'D430042O09Rik', 'Areg', 'Vmn2r90', 'Ier3', 'Cd14', 'Marcks', 'Angptl2', 'Rel', 'Prkg2', 'Afp', 'Serpinb2', 'Adh7', 'Btg2', 'Bdh2', 'Gm18445', 'Sdc4', 'Tnfsf9', 'Tpbg', 'Prr14l', 'Il27', 'Tk1', 'Angpt2', 'Tmem171', 'Ccl2', 'Ccl3', 'Ccl4', 'Sqstm1', 'Cd83', 'Slc7a11', 'Srl', 'Oasl1', 'Hsp90aa1', 'Slc9b2', 'Pde4b', 'Rasgrp3', 'Calcrl', 'Fosb', 'Egr1', 'Stx11', 'Colec12', 'Gmnn', 'Gpr84', 'Cxcl1', 'N28178', 'Cxcl2', 'Sod2', 'Zdhhc14', 'Mt2', 'Ttc7', 'Srgap3', 'Serpinb8', 'Srxn1', 'Phlda1', 'Apol9a', 'Bcl2a1d', 'Traf1', 'Serpinb9e', 'Pim1', 'Gm2004', 'Il1f9', 'Prdx1', 'Ccrl2', 'Slc1a2', 'Gem', 'Procr', 'Hmgcs1', 'Gm8325', 'AI607873', 'Bcl2l11', 'Col7a1', 'Irg1', 'Saa3', 'Gm28677', 'Tnf', 'Hdc', 'Arntl2', 'Sla', 'Il6', 'Dlg3', 'Klhl23', 'F8', 'Gpr183', 'Shc4', 'Fabp4', 'Atf3', 'Bco2', 'Ccl20', 'Ifit2', 'Errfi1', 'Lif', 'Kbtbd11', 'Sat1', 'Plaur', 'Bahcc1', 'Hmox1', 'Ak5', 'Id1', 'Mef2c', 'Icam1', 'Pcdh15', 'Slfn2', 'A930018M24Rik', 'Map1b', 'Lilrb4', 'Clec4e', 'Nfkbia', 'Csf3', 'Klhl6', 'Akr1b8', 'Emp1', 'Srgn', 'Zc3h12c', 'Prl2c2', 'Cela1', 'Slpi', 'Rasgef1b', 'Rnase4', 'Il23a', 'Mmp13', 'Plk3', 'Plk2', 'Rassf4', 'Stap1', 'Cish', 'Kdm6b', 'Il1a', 'Il1b', 'Gp49a', 'Malt1', 'Nabp1', 'Kif14', 'Ttc41', 'Rab31', 'Gsta1', 'Ppp1r15a', 'Hcar2', 'Myo18b', 'Mtss1', 'Ccr3', 'C130050O18Rik', 'Ccl5']
# print len(all_genes)
"""
Analysis!

Cluster the time traces and then compare the gene expression for each cluster
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
# import h5py
import pandas as pd
import numpy as np
import matplotlib as mpl
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2 import robjects as ro
import scipy.stats as stats

from dba import dba
from dba import align_to

from rpy2.robjects.vectors import DataFrame as RDataFrame
from rpy2 import rinterface
from rpy2.robjects import conversion

@conversion.py2ro.register(pd.DataFrame)
def py2ro_pandasdataframe(obj):
    ri_dataf = conversion.py2ri(obj)
    # cast down to an R list (goes through a different code path
    # in the DataFrame constructor, avoiding `str(k)`) 
    ri_list = rinterface.SexpVector(ri_dataf)
    return RDataFrame(ri_list)

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = .1
# mpl.style.use('ggplot')

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')
scde = importr("scde")

# Load data sets in R
r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")
r("""load("/scratch/PI/mcovert/dvanva/sequencing/counts_data.RData")""")
r("o.fpm = scde.expression.magnitude(o.ifm, counts = counts_data_int)")
list_of_all_genes = list(r("rownames(counts_data_int)"))
# list_of_genes = ['Wrb', 'Cxcl3', 'Tnfaip3', 'Tnfaip2', 'Blvrb', 'Nfatc1', 'Nfkbiz', 'Phf11d', 'Fas', 'Gpx1', 'Ccnd1', 'Rnf213', 'Marcksl1', 'Abcg1', 'Gm8818', 'Nlrp3', 'Rrm2', 'Mbnl1', 'Ptpn14', 'Odc1', 'Ptgs2', 'Ddit3', 'Cxcl10', 'Ly86', 'Ier3', 'Cd14', 'Rel', 'Prkg2', 'Afp', 'Btg2', 'Gm18445', 'Sdc4', 'Tnfsf9', 'Prr14l', 'Il27', 'Tk1', 'Angpt2', 'Tmem171', 'Ccl3', 'Ccl4', 'Sqstm1', 'Cd83', 'Slc7a11', 'Oasl1', 'Hsp90aa1', 'Pde4b', 'Rasgrp3', 'Calcrl', 'Egr1', 'Stx11', 'Colec12', 'Gmnn', 'Gpr84', 'Cxcl2', 'Sod2', 'Mt2', 'Serpinb8', 'Srxn1', 'Phlda1', 'Bcl2a1d', 'Traf1', 'Pim1', 'Il1f9', 'Prdx1', 'Procr', 'Hmgcs1', 'AI607873', 'Bcl2l11', 'Irg1', 'Saa3', 'Tnf', 'Hdc', 'Atf3', 'Errfi1', 'Lif', 'Sat1', 'Plaur', 'Hmox1', 'Id1', 'Mef2c', 'Icam1', 'Slfn2', 'Map1b', 'Lilrb4', 'Clec4e', 'Nfkbia', 'Csf3', 'Akr1b8', 'Emp1', 'Srgn', 'Zc3h12c', 'Slpi', 'Rasgef1b', 'Plk3', 'Plk2', 'Rassf4', 'Stap1', 'Kdm6b', 'Il1b', 'Gp49a', 'Malt1', 'Nabp1', 'Kif14', 'Rab31', 'Ppp1r15a', 'Mtss1', 'Ccl5']
# list_of_genes = ['B230307C23Rik', 'Wrb', 'Zfp677', 'Sele', 'Uox', 'Ptpn18', 'Gpr21', 'Fdft1', 'Akr1b8', 'Tmpo', 'Tfec', 'Zfhx4', 'Nova1', 'H2-M2', 'Ptger2', 'Gpx1', 'Il1f9', 'Rgs1', 'Gm14440', 'Gm14446', 'Il1f6', 'Rnf213', 'Kctd12', 'Mapk12', 'Glipr1', 'Idi1', 'Pde6b', 'Gm15459', 'Rrm2', 'Gdap1l1', 'Tnfaip2', 'Odc1', 'Fv1', 'Lrguk', 'Sytl3', 'Cebpd', 'Slc25a2', 'Ptprg', 'Gramd4', 'Tmem171', 'Maff', 'Ly86', 'D430042O09Rik', 'Areg', 'Ccnb1', 'Pim1', 'Ccdc17', 'Angptl2', 'Rel', 'Sprr2d', 'Cish', 'Gprc5a', 'Bcl6b', 'Zfp109', 'Mrc2', 'Btg2', 'Ryr3', 'Sdc4', 'Fabp4', 'Il27', 'Ifi44', 'Cd83', 'Gng8', 'Ccl2', 'Ccl3', 'Olfr810', 'Ccl4', 'Ccl5', 'Icos', 'Srl', 'Hist1h2ap', 'Hsp90aa1', 'Elmsan1', 'Sprr1a', 'Pde4b', 'Axl', 'Rasgrp3', 'Oscp1', 'AI429214', 'Calcrl', 'Lipn', 'Stx11', 'Colec12', 'Fas', 'Tnfsf4', 'Kifc3', 'Nkain1', 'Cxcl1', 'N28178', 'Cxcl2', 'Mt2', 'Slc17a6', 'Srgap3', 'P2rx3', 'Slc9b2', 'Dgkg', 'Olfr920', 'Sqstm1', 'Fut10', 'Gpr84', 'Pla2r1', 'Mturn', 'Hspa8', 'Kif18a', 'Zfp36', 'Plau', 'Ier3', 'Btnl9', 'Cd247', 'Gpr132', 'Ube2c', 'Gem', 'Dyrk4', 'Gins2', 'Hmgcs1', 'Gm8325', 'E2f3', 'Lcn2', 'Tnfsf18', 'Id1', 'F3', 'Lipg', 'Timp1', 'F8', 'Plek', 'Efna5', 'Gpr183', 'Fbxo15', 'Clcf1', 'Bco2', 'B4galt1', 'Ifit2', 'Errfi1', 'Lif', 'Kbtbd11', 'Sat1', 'Plaur', 'Cacna1b', 'Ak5', 'Il6', 'Il7', 'Gm6161', 'Icam1', 'Pcdh15', 'Slfn2', 'Map1a', 'Map1b', 'Lilrb4', 'Nlrp3', 'Gpr4', 'Clec4e', 'Csf1', 'Csf3', 'Klhl6', 'Csf3r', 'Emp1', 'Zc3h12c', 'Mmp3', 'Ccr1', 'Prob1', 'Rnase4', 'Mmp12', 'Mmp13', 'Wnk1', 'Rassf4', 'Wnk4', 'Ifi44l', 'Kdm6b', 'Pmaip1', 'Il1a', 'Il1b', 'Gm10767', 'Serpinb3b', 'Pla2g4e', 'Nr6a1', 'Kif14', 'Rab31', 'Il1rn', 'Fam214a', 'Birc5', 'Hcar2', 'Gm9919', 'Wdr25', 'Psd4', 'Gm597', 'Snx30', 'Ccr3', 'Ttc41', 'Rsad2', 'Ubb', 'Tm6sf2', 'Nfkbia', 'Gfi1b', 'Jam3', 'Usp2', 'Cxcl3', 'Dnmt3l', 'Ptpn14', 'Gm9687', 'Gm6377', 'Sprr2k', 'Mbnl1', 'Nfkb1', 'Nfkb2', 'Tcn2', 'Rsl1', 'Phf11d', 'Phf11a', 'Zbp1', 'Wdr95', 'Ccnd1', 'Flrt3', 'Gal', 'Spp1', 'Sirpb1a', 'Tslp', 'Atp6v0d2', 'Cytip', 'Car15', 'Atf3', 'Marcksl1', 'Abcg1', 'Gm8818', 'Ehd1', 'Nr4a3', 'Ubc', 'Tnfsf9', 'Zfp362', 'Cdk5r1', 'Myh14', 'Blvrb', 'Aars2', 'Ptgs2', 'Ddit3', 'Inpp5d', 'Fam189b', 'Slc6a9', 'Cfb', 'Col23a1', 'Plekho2', 'Cxcl10', 'Il4i1', 'Cela1', 'Accs', 'Gm13051', 'Kif7', 'Ifnb1', 'Npnt', 'Cd14', 'Marcks', 'Prkg2', 'Afp', 'Mipol1', 'Cdc6', 'Pdk1', 'Srxn1', 'Adh7', 'Serpinb2', 'Bdh2', 'Gm29633', 'Gm18445', 'Pglyrp3', 'Tpbg', 'Mamld1', 'Prr14l', 'Tk1', 'Angpt2', 'Ipcef1', 'Fsip2', 'Cysltr1', 'Shc4', 'Prl2c2', 'Slc3a2', 'Zfp708', '1700017B05Rik', 'Slc7a11', 'Oasl1', 'Per2', 'Tnfrsf11a', 'Ung', 'Bmf', 'Tmem132a', 'Pdcd1', 'Lgals3bp', 'Fosb', 'Egr1', 'Dclk1', 'Unc80', 'Gmnn', 'Ncmap', 'Txnrd1', 'Gm14419', 'Nupr1', 'Sod2', 'Tnfaip8l2', 'Zdhhc14', 'Gbp2', 'Epn1', 'Serpinb8', 'Zbtb42', 'Jun', 'Lpar1', 'Phlda1', 'Apol9b', 'Apol9a', 'Bcl2a1d', 'Traf1', 'Serpinb9b', 'Serpinb9e', 'Serpinb9g', 'Atp2a3', 'Lrp2bp', 'Vmn2r90', 'Gm2004', 'Prdx1', 'Ccrl2', 'Srgn', 'Slc1a2', 'Gm14401', 'Dusp16', 'Procr', 'Nfatc1', 'AI607873', 'Bcl2l11', 'Col7a1', 'St3gal4', 'Saa3', 'Slc9a9', 'Gm28677', 'Tnf', 'Hdc', 'Nabp1', 'Arntl2', 'P2ry6', 'Sla', 'Sema4g', 'Dlg3', 'Klhl23', 'Endod1', 'Tef', 'Zfp526', 'Slc4a5', 'Hmgb3', 'Slc4a7', 'Pdzk1ip1', 'Ccl20', 'Tnfaip3', 'Ccl22', 'Numbl', 'Ttc7', 'Siglech', 'Ankrd42', 'Erp29', 'Bahcc1', 'Siglecg', 'Hmox1', 'Id3', 'Gsta1', 'Mef2c', 'Mefv', 'Esd', 'Dusp5', 'A930018M24Rik', 'Dscaml1', 'Galnt3', 'Sqrdl', 'Cbx3', 'St6gal1', 'Hmga2', 'Nfkbiz', 'Myrfl', 'Olfr99', 'Edn1', 'Gadd45b', 'Slpi', 'Rasgef1b', 'Wscd2', 'Il23a', 'Plk1', 'Cdkl4', 'Plk3', 'Plk2', '2010111I01Rik', 'Gm4832', 'Stap1', 'Fam198b', 'Exoc3l2', 'Eps8', 'Insig1', 'Tnfrsf1b', 'Gp49a', 'Malt1', 'Irg1', 'Fbln2', 'Ccnb2', 'Adamts7', 'Cd40', 'Ppp1r15a', 'Myo18b', 'Mtss1', 'Pik3r3', 'C130050O18Rik', 'Tnfsf13b']
list_of_genes = ['3110043O21Rik', 'A430078G23Rik', 'AC102815.1', 'Aars', 'Acat2', 'Acot9', 'Actg1', 'Adh7', 'Adora2b', 'Adrm1', 'Afp', 'Akna', 'Akr1a1', 'Akr1b3', 'Akr1b8', 'Alas1', 'Ampd3', 'Anxa5', 'Anxa7', 'Areg', 'Arg2', 'Arih1', 'Arl14ep', 'Arl4a', 'Arl5b', 'Arl5c', 'Asns', 'Atf3', 'Atf4', 'Atp6v0b', 'BC028528', 'Bcat1', 'Bcl10', 'Bcl2a1b', 'Bcl2a1d', 'Bcl2l11', 'Bcl6b', 'Birc3', 'Blcap', 'Blvrb', 'Brd2', 'Brd4', 'Btg2', 'Bud31', 'C130026I21Rik', 'C920009B18Rik', 'Calcrl', 'Carhsp1', 'Cat', 'Ccdc8', 'Ccl2', 'Ccl20', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl9', 'Ccr3', 'Ccrl2', 'Cd14', 'Cd247', 'Cd40', 'Cd74', 'Cd83', 'Cdc42ep2', 'Cdc42ep4', 'Cdk18', 'Cdk5r1', 'Cdkn1a', 'Cebpd', 'Cfap61', 'Cflar', 'Chac1', 'Chordc1', 'Ciart', 'Cish', 'Cited2', 'Clec4d', 'Clec4e', 'Clic4', 'Cndp2', 'Creg1', 'Csf1', 'Csf3', 'Csrnp1', 'Cth', 'Ctsd', 'Cxcl10', 'Cxcl16', 'Cxcl2', 'Cxcl3', 'Cyb5a', 'Cybb', 'Cyp51', 'Dcbld2', 'Dcstamp', 'Ddit3', 'Ddx28', 'Dhrs3', 'Dnaja1', 'Dnajb1', 'Dnajb2', 'Dnajb4', 'Dnajb9', 'Dnmt3l', 'Dpep2', 'Dpp7', 'Dtx4', 'Dusp1', 'Dusp16', 'Dusp2', 'Dusp5', 'Dusp8', 'Dynll1', 'Ebi3', 'Eef1e1', 'Egr1', 'Egr2', 'Ehd1', 'Eif1', 'Eif3c', 'Eif4ebp1', 'Eif5', 'Ell2', 'Ercc1', 'Errfi1', 'Esd', 'Ets2', 'Evi2a', 'Exoc3l4', 'Fabp4', 'Fam53c', 'Fas', 'Fbxo15', 'Fbxo33', 'Fdft1', 'Fdps', 'Fgr', 'Flrt3', 'Fosl1', 'Fth1', 'Ftl1', 'Gabarap', 'Gabarapl1', 'Gadd45b', 'Galnt3', 'Gars', 'Gas2l3', 'Gbp2', 'Gch1', 'Gclm', 'Gdpd1', 'Gem', 'Ghitm', 'Glipr1', 'Glrx', 'Gm10116', 'Gm10718', 'Gm10800', 'Gm11425', 'Gm11427', 'Gm13502', 'Gm15459', 'Gm1821', 'Gm18445', 'Gm29013', 'Gm4832', 'Gm5526', 'Gm8818', 'Gm8995', 'Gm9687', 'Gm9938', 'Gmfb', 'Gp49a', 'Gpatch11', 'Gpr132', 'Gpr84', 'Gramd1b', 'Gsap', 'Gspt1', 'Gsta1', 'Gtf2b', 'H2-Q4', 'H2-Q7', 'H60b', 'Hax1', 'Hcar2', 'Hdc', 'Herpud1', 'Hilpda', 'Hmgcr', 'Hmgcs1', 'Hmox1', 'Hsd17b7', 'Hsp90aa1', 'Hsp90ab1', 'Hspa1b', 'Hspa8', 'Hsph1', 'Htatip2', 'Icam1', 'Id1', 'Id3', 'Idi1', 'Ier3', 'Ier5', 'Ifih1', 'Ifrd1', 'Igsf3', 'Igsf6', 'Il18rap', 'Il1a', 'Il1b', 'Il1f6', 'Il1f9', 'Il1rn', 'Il27', 'Il2rg', 'Il4i1', 'Il6', 'Il7', 'Impact', 'Insig1', 'Ipo13', 'Irf1', 'Irg1', 'Isg15', 'Itga5', 'Jag1', 'Jak2', 'Jun', 'Junb', 'Kansl1l', 'Kcnj2', 'Kctd12', 'Kdm6b', 'Kif3c', 'Klf6', 'Kpna4', 'Lamc2', 'Layn', 'Lcn2', 'Ldlr', 'Lif', 'Lilrb4', 'Litaf', 'Lpar1', 'Lpin1', 'Lss', 'Mad2l1bp', 'Maff', 'Malt1', 'Map1b', 'Map1lc3b', 'Mapk6', 'Marcks', 'Marcksl1', 'Mcl1', 'Med10', 'Mettl4', 'Mitf', 'Mllt11', 'Mmp13', 'Mndal', 'Mrpl52', 'Msantd3', 'Msmo1', 'Mt1', 'Mt2', 'Mthfd2', 'Mtmr7', 'Mvd', 'Nabp1', 'Narf', 'Nars', 'Nfe2l1', 'Nfil3', 'Nfkb1', 'Nfkb2', 'Nfkbia', 'Nfkbib', 'Nfkbid', 'Nfkbie', 'Nfkbiz', 'Ninj1', 'Nlrp3', 'Npnt', 'Nqo1', 'Nr4a3', 'Nsdhl', 'Nub1', 'Nupr1', 'Oasl1', 'Odc1', 'Oser1', 'Osgin1', 'Osgin2', 'P2rx4', 'Parp3', 'Pcdh15', 'Pdcd1', 'Pde4b', 'Pdgfb', 'Pdlim7', 'Pdp2', 'Pea15a', 'Pef1', 'Peli1', 'Pgd', 'Phlda1', 'Pik3r5', 'Pim1', 'Pinx1', 'Pla2g4e', 'Plagl2', 'Plaur', 'Plek', 'Plekhm3', 'Plekho2', 'Plin2', 'Plk2', 'Plk3', 'Plscr1', 'Pmaip1', 'Pmvk', 'Ppp1r11', 'Ppp1r15a', 'Ppp4r2', 'Prdx1', 'Prdx5', 'Prdx6', 'Prkg2', 'Procr', 'Prr13', 'Psmb6', 'Psmc2', 'Psmc6', 'Psmd10', 'Psmd11', 'Psmd7', 'Psmd8', 'Psph', 'Ptger2', 'Ptgr1', 'Ptgs2', 'Ptpn14', 'Ptpre', 'Pvr', 'RP23-172K2.12', 'RP23-57A17.2', 'Rab11fip1', 'Rab32', 'Rab8b', 'Rap1b', 'Rap2b', 'Rasgef1b', 'Rassf4', 'Rc3h1', 'Rdh11', 'Rel', 'Relb', 'Rffl', 'Rgs1', 'Rgs16', 'Rhbdf2', 'Rhob', 'Rhoc', 'Rhof', 'Riok3', 'Rnf19b', 'Rrs1', 'Saa3', 'Sars', 'Sat1', 'Sc5d', 'Sdc4', 'Sde2', 'Sele', 'Selk', 'Sept11', 'Serpinb1b', 'Serpinb2', 'Serpinb9b', 'Serpinb9e', 'Serpine1', 'Sertad1', 'Sertad2', 'Sgk1', 'Skil', 'Slc11a1', 'Slc15a3', 'Slc17a6', 'Slc1a4', 'Slc25a33', 'Slc25a37', 'Slc2a6', 'Slc39a1', 'Slc3a2', 'Slc48a1', 'Slc4a7', 'Slc6a9', 'Slc7a11', 'Slc7a6os', 'Slc9a4', 'Slfn10-ps', 'Slfn2', 'Slpi', 'Snx10', 'Socs3', 'Socs5', 'Sod2', 'Spp1', 'Spryd7', 'Spty2d1', 'Sqrdl', 'Sqstm1', 'Srfbp1', 'Srgn', 'Srxn1', 'Stap1', 'Stard4', 'Stip1', 'Stk40', 'Stx11', 'Stx6', 'Taf13', 'Taf7', 'Tank', 'Tars', 'Tbpl1', 'Tbrg1', 'Tfec', 'Tgif1', 'Tgoln1', 'Tiparp', 'Tlcd2', 'Tlr2', 'Tma16', 'Tmbim6', 'Tmem171', 'Tmem41b', 'Tmsb10', 'Tnf', 'Tnfaip2', 'Tnfaip3', 'Tnfaip8l2', 'Tnfrsf10b', 'Tnfrsf1b', 'Tnfsf4', 'Tnfsf9', 'Tnip1', 'Tnip3', 'Top1', 'Tox4', 'Tpbg', 'Tpm4', 'Traf1', 'Treml4', 'Trex1', 'Trib1', 'Trim13', 'Trpv2', 'Tslp', 'Tubb2a', 'Txndc9', 'Txnrd1', 'Uba6', 'Ubb', 'Ubc', 'Ube2h', 'Ubxn4', 'Unc80', 'Vegfa', 'Vmn2r90', 'Wsb1', 'Ypel5', 'Zbtb32', 'Zc3h12a', 'Zc3h12c', 'Zc3hav1', 'Zfand2a', 'Zfand5', 'Zfp330', 'Zfp36', 'Zfp617', 'Zfp622', 'Znrf1', 'Zswim4', 'Zwint', 'Zyx']

# r("write.csv(exp(o.fpm), file = '/scratch/PI/mcovert/dvanva/sequencing/fpms.csv')")

fpm = pd.read_csv('/scratch/PI/mcovert/dvanva/sequencing/fpms.csv', index_col = 0)
logfpm = np.log2(fpm + 1)

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'
all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))

all_cells = []
longest_time = 0

time_point = 300
cell_names = []
for cell in all_cells_total:
	if cell.time_point == 75 or cell.time_point == 150 or cell.time_point == 300 and cell.condition == 'Stim':
		longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
		all_cells += [cell]
		cell_names += [cell.id]

reduced_logfpm = logfpm.loc[:,cell_names]

gene1_array = reduced_logfpm.loc['Ccl3',:].dropna()
gene2_array = reduced_logfpm.loc['Ccl4',:].dropna()

pear = scipy.stats.pearsonr(gene1_array, gene2_array)[0]
print pear
"""
Plot dependence of average fpm on dropout rate
"""

avg_fpm = []
dropout_rate = []
gene_names = []
coeff_var = []

counter = 0
total_cells = np.float(len(all_cells))
# print total_cells
for gene in list_of_genes:
	fpms = reduced_logfpm.loc[gene,:]
	dropout_rate += [1-np.float(fpms.nonzero()[0].shape[0])/total_cells]
	avg_fpm +=[np.mean(fpms.iloc[fpms.nonzero()])]
	gene_names += [gene]
	coeff_var += [scipy.stats.variation(fpms.dropna())]

# print dropout_rate
# print avg_fpm

fig = plt.figure(figsize = (5,4))
ax = fig.add_subplot(111)
ax.scatter(avg_fpm, dropout_rate, color = 'b', s = .5, alpha = 1)
# ax.set_xscale('log')
ax.set_xlabel('Mean log_2(FPM+1) across all cells', fontsize = 12)
ax.set_ylabel('Dropout frequency', fontsize = 12)
ax.set_title('Dropout frequency vs expression', y= 1.05, fontsize = 14)
ax.set_ylim([-0.05,1.05])
ax.set_yticks([0,.2,.4,.6,.8,1])
ax.set_xlim([0,20])

for counter in xrange(len(gene_names)):
	ax.annotate(gene_names[counter], xy = (avg_fpm[counter], dropout_rate[counter]), textcoords = 'data', fontsize = 3)

fig.tight_layout()
plt.savefig("plots/trial_30_dropout_freq.pdf")

"""
Plot coefficient of variation
"""

plt.gcf()
fig = plt.figure(figsize = (5,4))
ax = fig.add_subplot(111)

ind = np.arange(len(coeff_var))
width = 0.01
rects = ax.scatter(avg_fpm, coeff_var, color = 'b', s = .5)
# ax.set_xticks(ind + width)
# ax.set_xticklabels(gene_names, fontsize = 3, rotation = 90)
ax.set_ylabel('Coefficient of variation of log2(FPM+1)', fontsize = 12)
ax.set_xlabel('Mean log2(FPM+1) across all cells', fontsize = 12)
ax.set_xlim([0,20])

# ax.set_ylim([-0.05,1.05])
# ax.set_yticks([0,.2,.4,.6,.8,1])
for counter in xrange(len(gene_names)):
	if gene_names[counter] == 'Ccl5' or gene_names[counter] == 'Cxcl10':
		ax.annotate(gene_names[counter], xy = (avg_fpm[counter], coeff_var[counter]), textcoords = 'data', fontsize = 3, color = 'r')
	else:
		ax.annotate(gene_names[counter], xy = (avg_fpm[counter], coeff_var[counter]), textcoords = 'data', fontsize = 3)

fig.tight_layout()
plt.savefig("plots/trial_30_coeff_var.pdf")

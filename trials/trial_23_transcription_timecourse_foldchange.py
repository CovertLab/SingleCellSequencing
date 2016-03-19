

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

from dba import dba
from dba import align_to

from rpy2.robjects.vectors import DataFrame as RDataFrame
from rpy2 import rinterface
from rpy2.robjects import conversion
import seaborn as sns

@conversion.py2ro.register(pd.DataFrame)
def py2ro_pandasdataframe(obj):
    ri_dataf = conversion.py2ri(obj)
    # cast down to an R list (goes through a different code path
    # in the DataFrame constructor, avoiding `str(k)`) 
    ri_list = rinterface.SexpVector(ri_dataf)
    return RDataFrame(ri_list)

def zero_crossing(data, offset):
	new_data = data - offset
	zc = np.where(np.diff(np.signbit(new_data)))[0]
	return zc[0]

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
# mpl.style.use('ggplot')

R = rpy2.robjects.r
DTW = importr('dtw')
DTWclust = importr('dtwclust')
scde = importr("scde")


# Load data sets in R
r("""load("/scratch/PI/mcovert/dvanva/sequencing/all_cells_scde_fit_linear.RData")""")
r("""load("/scratch/PI/mcovert/dvanva/sequencing/counts_data.RData")""")
list_of_all_genes = list(r("rownames(counts_data_int)"))

# Load pickle file with cell objects
direc = '/scratch/PI/mcovert/dvanva/sequencing/'
all_cell_file = 'all_cells_qc_complete.pkl'

all_cells_total = pickle.load(open(os.path.join(direc,all_cell_file)))


# Determine which genes to look at
# list_of_genes = ["Nfkbia", "Tnfaip3", "Nfkbie", "Nfkbiz", "Nlrp3", "Tnfaip2", "Atf3"]
# list_of_genes = ["Cxcl2", "Cxcl3", "Cxcl10", "Ccl3", "Ccl4", "Ccl5", "Ccl20", "Tnf", "Tnfsf9", "Il1a", "Il1b", "Il1f6", "Il6", "Il27", "Il1f9", "Lif",  "Csf3"]
# list_of_genes = ["Hmox1", "Prdx1", "Hdc", "Ptgs2", "Irg1"]
# list_of_genes = ["Plaur", "Sqstm1", "Clec4e", "Sdc4", "Procr", "Slpi", "Plk2", "Saa3", "Slc7a11", "Cish", "Gp49a", "Hcar2", "Gpr84", "Malt1"]
# list_of_genes = ['Dnmt3b', 'Tecta', 'Tm6sf2', 'Bricd5', 'Prdm12', 'Prdm13', 'Adora2a', 'Ccdc162', 'Gm5283', 'Gm11400', 'Olfr536', 'Gm13145', 'Gm13333', 'Zfp661', 'Angptl3', 'Sipa1l3', 'Scn1a', 'Sprr2d', 'Il17rc', 'Zglp1', 'Akr1cl', 'Map1a', 'Trim14', 'Adgrg6', 'Gm13991', 'Dhrs11', 'Gm21834', 'Iqca', 'Gm2007', 'Slc39a8', 'Gng7', 'AL663030.1', 'Nphp4', 'Nod1', 'Emc9', 'Akr1b7', 'Il33', 'Mmp14', 'Zfyve1', 'Cetn4', '2610305D13Rik', 'Mettl25', 'Ric8b', 'Mterf2', 'Zfp850', 'Clec4a4', 'Saa3', 'Hist1h4n', 'Gm11007', 'Cntrob', 'Atp7b', 'Mtl5', '1700061G19Rik', 'Coro2b', '1700030J22Rik', 'Gm8898', 'Tmem86b', 'Car9', 'Gm5157', 'Gm15539', 'Arhgef18', 'Slc13a3', 'Dclk1', 'Ager', 'Actr3b', 'Zfp41', 'Fzd8', '4930524J08Rik', 'Zic5', 'Trem1', 'Ppp1r32', 'Stk36', 'Gnao1', 'Tmem239', 'Polm', 'Fgf21', 'Gprasp2', 'Tesk1', 'Athl1', 'Kptn']
# list_of_genes = ['Mmp3', 'Ccl5', 'Gpr137c', 'Efna5', 'Tiam1', 'D2hgdh', 'Nod2', 'Gm14440', 'Pla2r1', 'Serpinb9g', 'Hic2', 'Cdkl4', 'Slc18b1', 'H2-M2', 'Klhdc1', 'Iqcb1', 'Sh3bp2', 'Ifit3', 'Cmpk2', 'Adamts10', 'Sirt5', 'Plekhg2', 'Cxcl10', 'Gm13051', 'Tppp3', 'Krt24', 'Lamb3', 'Serpind1', 'Pars2', 'Spopl', 'Rsad2', 'Tnfsf4', 'Gm12728', 'Siglece', '4930432K21Rik', 'Vmn1r32', 'Fbxw10', 'Ngb', 'Bdkrb1', 'B3galt2']
# list_of_genes = ["Saa3", "Il1f9", "Il6", "Csf3", "Gp49a"]
# list_of_genes = ["Cxcl2", "Cxcl10", "Ccl20", "Tnfsf9", "Il1a"]
# list_of_genes = ["Tnf", "Ccl4", "Il1b", "Cxcl3", "Il1f6", "Lif"]
# list_of_genes = ["Ccl3", "Ccl5"]
# list_of_genes = ["Saa3", "Il1f9", "Il6", "Csf3", "Gp49a", "Tm6sf2", "Il17rc", "Zfyve1", "Ric8b", "Sipa1l3", "Zfp850", "Nod1", "Trim14", "Trem1", "Athl1"]
# list_of_genes = ["Prdx1", "Irg1"]
# list_of_genes = ['Wrb', 'Cxcl3', 'Tnfaip3', 'Tnfaip2', 'Gm6377', 'Blvrb', 'Bcl6b', 'Gpr21', 'Rsl1', 'Nfatc1', 'Nfkbiz', 'Phf11d', 'Tfec', 'Phf11a', 'Fas', 'Gpx1', 'Ccnd1', 'Flrt3', 'Cd247', 'Tslp', 'Atp6v0d2', 'Il1f6', 'Rnf213', 'Marcksl1', 'Abcg1', 'Gm8818', 'Nlrp3', 'Pde6b', 'Rrm2', 'Mbnl1', 'Ptpn14', 'Odc1', 'Fv1', 'Ptgs2', 'Ddit3', 'Slc25a2', 'Cfb', 'Mmp3', 'Ptprg', 'Cxcl10', 'Ly86', 'D430042O09Rik', 'Areg', 'Vmn2r90', 'Ier3', 'Cd14', 'Marcks', 'Angptl2', 'Rel', 'Prkg2', 'Afp', 'Serpinb2', 'Adh7', 'Btg2', 'Bdh2', 'Gm18445', 'Sdc4', 'Tnfsf9', 'Tpbg', 'Prr14l', 'Il27', 'Tk1', 'Angpt2', 'Tmem171', 'Ccl2', 'Ccl3', 'Ccl4', 'Sqstm1', 'Cd83', 'Slc7a11', 'Srl', 'Oasl1', 'Hsp90aa1', 'Slc9b2', 'Pde4b', 'Rasgrp3', 'Calcrl', 'Fosb', 'Egr1', 'Stx11', 'Colec12', 'Gmnn', 'Gpr84', 'Cxcl1', 'N28178', 'Cxcl2', 'Sod2', 'Zdhhc14', 'Mt2', 'Ttc7', 'Srgap3', 'Serpinb8', 'Srxn1', 'Phlda1', 'Apol9a', 'Bcl2a1d', 'Traf1', 'Serpinb9e', 'Pim1', 'Gm2004', 'Il1f9', 'Prdx1', 'Ccrl2', 'Slc1a2', 'Gem', 'Procr', 'Hmgcs1', 'Gm8325', 'AI607873', 'Bcl2l11', 'Col7a1', 'Irg1', 'Saa3', 'Gm28677', 'Tnf', 'Hdc', 'Arntl2', 'Sla', 'Il6', 'Dlg3', 'Klhl23', 'F8', 'Gpr183', 'Shc4', 'Fabp4', 'Atf3', 'Bco2', 'Ccl20', 'Ifit2', 'Errfi1', 'Lif', 'Kbtbd11', 'Sat1', 'Plaur', 'Bahcc1', 'Hmox1', 'Ak5', 'Id1', 'Mef2c', 'Icam1', 'Pcdh15', 'Slfn2', 'A930018M24Rik', 'Map1b', 'Lilrb4', 'Clec4e', 'Nfkbia', 'Csf3', 'Klhl6', 'Akr1b8', 'Emp1', 'Srgn', 'Zc3h12c', 'Prl2c2', 'Cela1', 'Slpi', 'Rasgef1b', 'Rnase4', 'Il23a', 'Mmp13', 'Plk3', 'Plk2', 'Rassf4', 'Stap1', 'Cish', 'Kdm6b', 'Il1a', 'Il1b', 'Gp49a', 'Malt1', 'Nabp1', 'Kif14', 'Ttc41', 'Rab31', 'Gsta1', 'Ppp1r15a', 'Hcar2', 'Myo18b', 'Mtss1', 'Ccr3', 'C130050O18Rik', 'Ccl5']
# list_of_genes = ['Wrb', 'Cxcl3', 'Tnfaip3', 'Tnfaip2', 'Blvrb', 'Nfatc1', 'Nfkbiz', 'Phf11d', 'Fas', 'Gpx1', 'Ccnd1', 'Rnf213', 'Marcksl1', 'Abcg1', 'Gm8818', 'Nlrp3', 'Rrm2', 'Mbnl1', 'Ptpn14', 'Odc1', 'Ptgs2', 'Ddit3', 'Cxcl10', 'Ly86', 'Ier3', 'Cd14', 'Rel', 'Prkg2', 'Afp', 'Btg2', 'Gm18445', 'Sdc4', 'Tnfsf9', 'Prr14l', 'Il27', 'Tk1', 'Angpt2', 'Tmem171', 'Ccl3', 'Ccl4', 'Sqstm1', 'Cd83', 'Slc7a11', 'Oasl1', 'Hsp90aa1', 'Pde4b', 'Rasgrp3', 'Calcrl', 'Egr1', 'Stx11', 'Colec12', 'Gmnn', 'Gpr84', 'Cxcl2', 'Sod2', 'Mt2', 'Serpinb8', 'Srxn1', 'Phlda1', 'Bcl2a1d', 'Traf1', 'Pim1', 'Il1f9', 'Prdx1', 'Procr', 'Hmgcs1', 'AI607873', 'Bcl2l11', 'Irg1', 'Saa3', 'Tnf', 'Hdc', 'Atf3', 'Errfi1', 'Lif', 'Sat1', 'Plaur', 'Hmox1', 'Id1', 'Mef2c', 'Icam1', 'Slfn2', 'Map1b', 'Lilrb4', 'Clec4e', 'Nfkbia', 'Csf3', 'Akr1b8', 'Emp1', 'Srgn', 'Zc3h12c', 'Slpi', 'Rasgef1b', 'Plk3', 'Plk2', 'Rassf4', 'Stap1', 'Kdm6b', 'Il1b', 'Gp49a', 'Malt1', 'Nabp1', 'Kif14', 'Rab31', 'Ppp1r15a', 'Mtss1', 'Ccl5']
# add_genes = ['Nfkbie', 'Il6', 'Ccl20', 'Il1a', 'Il1f6', 'Zfp850']
# list_of_genes = list(set(list_of_genes)|set(add_genes))


# list_of_genes = ['Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl6', 'Ccl7', 'Ccl9', 'Ccl11', 'Ccl17', 'Ccl20', 'Ccl22', 'Ccl24', 'Ccl25', 'Ccl28']
# list_of_genes = ['Tlr1', 'Tlr2', 'Tlr4', 'Tlr6', 'Cd14', 'Ly96'] #All TLR genes
# list_of_genes = ['Rela', 'Zfp36', 'Ifnb1']
# list_of_genes = ['Egr1', 'Fos', 'Nr4a1', 'Btg2', 'Zfp36', 'Fabp4', 'C5ar1', 'Ier3', 'Gdf15', 'Maff', 'Cd83', 'Egr2', 'Gpr84', 'Dusp2', 'Sdc4', 'Errfi1', 'Dusp5', 'Zc3h12c', 'Plek', 'Kdm6b', 'Nlrp3', 'Ifrd1', 'Dusp4', 'Rasgef1b','Jag1', 'Mapkapk2', 'Rapgef2', 'Sqstm1', 'Trib1', 'Rgs1', 'Cd14', 'Cxcl1', 'Il1b', 'Cxcl2', 'Ptgs2', 'Tnfsf9', 'Nfkbiz'] #MAPK Dependent
# list_of_genes = ['Ccl5', 'Cxcl10', 'Gbp5', 'Irg1', 'Isg15', 'Ifnb1', 'Ifih1', 'Peli1', 'Slc25a37'] #IRF3 Dependent
# list_of_genes = ['Il1rn', 'Casp4', 'Ccl2', 'Ccl7', 'Stx11', 'Ripk2', 'Gm6377', 'Cybb', 'Prdm1', 'Cd44', 'Ccl9', 'St3gal1', 'Il1a', 'Abtb2', 'Slc7a2', 'Malt1', 'Zhx2', 'Gpr85', 'Arhgef3', 'Fnbp1l', 'Hivep2', 'Slc44a1', 'Fchsd2', 'Aoah', 'Parp9', 'Ube2e2', 'Rnd3'] #TRIF Dependent

# list_of_genes = ['Tnfsf15','Tabpbp','Dtx3I','Gm12216','Parp14','Il15','Trim21','Herc6','Ccnd2','AA467197','Il15ra','Gbp7','Ms4a6c','Mndal','Ifi203','Pnp','Ifi205','Nod1','Cd274','Gnb4','Oasl1','Nt5c3','Ppm1k','Cxcl11','Ifi204','Slfn8','Ifit1','Batf2','Cxcl9','Rsad2','Ifi35','Oasl2','Isg20','Ddx60','Cmpk2','Trim30a','Tor3a','Zbp1','Mx1','Mx2','Daxx','Tnfsf10','Ifit3','Tap1','Parp12','Pyhin1','Igtp','Ifi47','Ifit2','Olfr56','Gm12250','Slfn5','Usp18','Irgm2','Trim30a','Bst2','Irf7','D14Ertd668e','Setdb2','Irgm1','Xaf1','Stat1','Ddx58','Oas1b','Oas1g'] # INFAR dependent
# list_of_genes = ['IL4i1','Saa3','D8Ertd82e','Adora2a','Stat5a','Jdp2','Il12b','Slc7a11','Lass6','Cav1','Jak2','Pstpip2','Fmnl2','Mmp13','Nupr1','Cish','Slco3a1','Snn','C3','Lif','Il27','Psme2','Clec2d','Il6','Nos2','Il18','Gbp2','Gbp3','Cd69'] #IFNAR independent
# list_of_genes=['Il1m','Ripk2','Gm6377','Il1a','Ccl9','Aoah','Rnd3','Errfi1','Il1b','Ptgs2','Med21','Csf1','Traf1','Tnip3','Ebi3'] #Myd88 dependent
# list_of_genes = ['Clec4e','Pde4b','Sod2','Socs3','Srgn','Gem','Marcksl1','Tnf','Ehd1','Serpine1','Csmp1','Dusp2','Il1b','Tnfsf9','Nlrp3','Plek','Zc3h12c','Maff','Mapkapk2','Jag1','Errfi1','Rasgef1b','Sqstm1'] # Nfkb enhancer
# list_of_genes = ['Tlr2','Nfkbib','Sdc4','Relb','Nfkbid','Bcl2l11','Ebi3','Cxcl2','Cxcl1','Cxcl10','Ifnb1','Tnfaip3','Casp4','Nfkb1','Nfkbia','Nfkbiz','Fchsd2','Gbp5','Ptgs2','Nfkbie','Nfkb2','Csf1','Stx11','Rapgef2','Ier3','Irg1','Kdm6b','Tnip3','Cd40','Icam1','Traf1','Irf1','Ccl5','Gpr84','Btg2','Ccrl2','Cd83'] # Strong RelA motif and chip peak
# list_of_genes = ['Serpine1','Agrn','Phldb1','Ccl9','Hivep2','Src','Sod2'] # rela motif no chip peak
# list_of_genes = ['Nfkbid','Nfkbia','Tlr2','Nfkbib','Nfkbie','Relb','Nfkb2','Tnfaip3','Tnip3','Traf1','Cd40','Nfkb1','Nfkbiz','Btg2','Ier3','Cd83','Gpr84','Sdc4','Kdm6b','Rapgef2','Cxcl2','Cxcl1','Ptgs2','Irf1','Ebi3','Bcl2l11','Ccrl2','Icam1','Csf1','Stx11','Casp4','Fchsd2','Ccl5','Cxcl10','Irg1','Gbp5','Ifnb1','Clec4e','Pde4b','Sod2','Socs3','Srgn','Gem','Marcksl1','Tnf','Ehd1','Serpine1','Csrnp1','Dusp2','Il1b','Tnfsf9','Nlrp3','Plek','Zc3h12c','Egr2','Maff','Mapkapk2','Jag1','Errfi1','Dusp5','Rasgef1b','Sqstm1'] # NFkb targets
# list_of_genes = ['Ccl5','Cxcl10','Gbp5','Irg1','Isg15','Ifnb1','Il1m','Casp4','Tcfec','Ccl2','Ccl7','Gm6377','Cybb','Cd44','Ccl9','Il1a','Arhgef3','Aoah','Fabp4','C5ar1','Gpr84','Plek','Kdm6b','Nlrp3','Rapgef2','Rgs1','Cd14','Il1b','Cxcl2','Ccl4','Ccl3','Slfn2','Arl5c','Ebi3','Cxcl16','Ccrl2','Srgn','Tnf','Clec4e','Pde4b','Tnip3','Serpine1','Cd40','Traf1','BC006779','Atp2b4'] #CpG low
# list_of_genes = ['Ifih1','Peli1','Slc25a37','Stx11','Ripk2','Prdm1','St3gal1','Abtb2','Slc7a2','Malt1','Zhx2','Gpr85','Fnbp1l','Hivep2','Slc44a1','Fchsd2','Parp9','Ube2e2','Rnd3','Fos1','Nr4a1','Dusp1','Btg2','Zfp36','Ier3','Gdf15','Maff','Cd83','Egr2','Dusp2','Sdc4','Errfi1','Dusp5','Zc3h12c','Ifrd1','Dusp4','Rasgef1b','Jag1','Mapkapk2','Sqstm1','Trib1','Cxcl1','Ptgs2','Tnfsf9','Nfkbiz','Ppp1r15a','Ier2','Nfkbid','Irf1','Nfkbia','Tlr2','Med21','Nfkbib','Rnf19b','Ccrn4l','Nfkbie','Frmd6','Icosl','Tnfaip2','Orai2','Slc2a6','Bcl2l11','Rcan1','Relb','Clic4','Nfkb2','Cflar','Gem','Socs3','Tnfaip3','Marcksl1','Csrnp1','Sod2','Ehd1','Csf1','Icam1','Src','Phldb1','Nod2','Rab11fip1','Foxp4','Flnb','Nfkb1','Agrn','Col18a1'] #CpG genes
# list_of_genes = ['Serpine1','Agrn','Phldb1','Ccl9','Hivep2','Src','Sod2'] # strong rela motif no chip peak
# list_of_genes = ['Zhx2','Tcfec','Cxcl16','Dusp5','Ier2','Cd44','Rasgef1b','Il1m','Il1b','Slc2a6','Srgn','Ccl3','Tnfsf9','Nr4a1','Cflar','Gm6377','Gem','Ccl4','Cybb','Slfn2','Malt1','Tnfsf9','C5ar1','Plek','Clec4e','Aoah'] # weak rela motif and chip peak
# list_of_genes = ['Btg2','Ier3','Cd83','Gpr84','Sdc4','Kdm6b','Rapgef2','Cxcl2','Cxl1','Ptgs2'] # Nfkb binding site and mapk dependent
# list_of_genes = ['Ccl5','Cxcl10','Irg1','Gbp5','Ifnb1'] # Nfkb site IRF3 dependent
# list_of_genes = ['Ccl5','Cxcl10','Gbp5','Irg1','Isg15','Ifnb1','Ifih1','Peli1','Slc25a37','Ccrn4l','Agrn','Cxcl1','C5ar1','Parp9','Slc7a2','Ccl9'] #IRF3 dependent all
# list_of_genes = ['Rnd3','Nr4a1','Egr1','Fos','Zfp36','Egr2','Dusp5'] #SRF
# list_of_genes = ['Tnfsf10', 'Ccl2', 'Casp4', 'Serpina3g', 'Ripk2']
# list_of_genes = ['3110043O21Rik', 'A430078G23Rik', 'AC102815.1', 'Aars', 'Acat2', 'Acot9', 'Actg1', 'Adh7', 'Adora2b', 'Adrm1', 'Afp', 'Akna', 'Akr1a1', 'Akr1b3', 'Akr1b8', 'Alas1', 'Ampd3', 'Anxa5', 'Anxa7', 'Areg', 'Arg2', 'Arih1', 'Arl14ep', 'Arl4a', 'Arl5b', 'Arl5c', 'Asns', 'Atf3', 'Atf4', 'Atp6v0b', 'BC028528', 'Bcat1', 'Bcl10', 'Bcl2a1b', 'Bcl2a1d', 'Bcl2l11', 'Bcl6b', 'Birc3', 'Blcap', 'Blvrb', 'Brd2', 'Brd4', 'Btg2', 'Bud31', 'C130026I21Rik', 'C920009B18Rik', 'Calcrl', 'Carhsp1', 'Cat', 'Ccdc8', 'Ccl2', 'Ccl20', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl9', 'Ccr3', 'Ccrl2', 'Cd14', 'Cd247', 'Cd40', 'Cd74', 'Cd83', 'Cdc42ep2', 'Cdc42ep4', 'Cdk18', 'Cdk5r1', 'Cdkn1a', 'Cebpd', 'Cfap61', 'Cflar', 'Chac1', 'Chordc1', 'Ciart', 'Cish', 'Cited2', 'Clec4d', 'Clec4e', 'Clic4', 'Cndp2', 'Creg1', 'Csf1', 'Csf3', 'Csrnp1', 'Cth', 'Ctsd', 'Cxcl10', 'Cxcl16', 'Cxcl2', 'Cxcl3', 'Cyb5a', 'Cybb', 'Cyp51', 'Dcbld2', 'Dcstamp', 'Ddit3', 'Ddx28', 'Dhrs3', 'Dnaja1', 'Dnajb1', 'Dnajb2', 'Dnajb4', 'Dnajb9', 'Dnmt3l', 'Dpep2', 'Dpp7', 'Dtx4', 'Dusp1', 'Dusp16', 'Dusp2', 'Dusp5', 'Dusp8', 'Dynll1', 'Ebi3', 'Eef1e1', 'Egr1', 'Egr2', 'Ehd1', 'Eif1', 'Eif3c', 'Eif4ebp1', 'Eif5', 'Ell2', 'Ercc1', 'Errfi1', 'Esd', 'Ets2', 'Evi2a', 'Exoc3l4', 'Fabp4', 'Fam53c', 'Fas', 'Fbxo15', 'Fbxo33', 'Fdft1', 'Fdps', 'Fgr', 'Flrt3', 'Fosl1', 'Fth1', 'Ftl1', 'Gabarap', 'Gabarapl1', 'Gadd45b', 'Galnt3', 'Gars', 'Gas2l3', 'Gbp2', 'Gch1', 'Gclm', 'Gdpd1', 'Gem', 'Ghitm', 'Glipr1', 'Glrx', 'Gm10116', 'Gm10718', 'Gm10800', 'Gm11425', 'Gm11427', 'Gm13502', 'Gm15459', 'Gm1821', 'Gm18445', 'Gm29013', 'Gm4832', 'Gm5526', 'Gm8818', 'Gm8995', 'Gm9687', 'Gm9938', 'Gmfb', 'Gp49a', 'Gpatch11', 'Gpr132', 'Gpr84', 'Gramd1b', 'Gsap', 'Gspt1', 'Gsta1', 'Gtf2b', 'H2-Q4', 'H2-Q7', 'H60b', 'Hax1', 'Hcar2', 'Hdc', 'Herpud1', 'Hilpda', 'Hmgcr', 'Hmgcs1', 'Hmox1', 'Hsd17b7', 'Hsp90aa1', 'Hsp90ab1', 'Hspa1b', 'Hspa8', 'Hsph1', 'Htatip2', 'Icam1', 'Id1', 'Id3', 'Idi1', 'Ier3', 'Ier5', 'Ifih1', 'Ifrd1', 'Igsf3', 'Igsf6', 'Il18rap', 'Il1a', 'Il1b', 'Il1f6', 'Il1f9', 'Il1rn', 'Il27', 'Il2rg', 'Il4i1', 'Il6', 'Il7', 'Impact', 'Insig1', 'Ipo13', 'Irf1', 'Irg1', 'Isg15', 'Itga5', 'Jag1', 'Jak2', 'Jun', 'Junb', 'Kansl1l', 'Kcnj2', 'Kctd12', 'Kdm6b', 'Kif3c', 'Klf6', 'Kpna4', 'Lamc2', 'Layn', 'Lcn2', 'Ldlr', 'Lif', 'Lilrb4', 'Litaf', 'Lpar1', 'Lpin1', 'Lss', 'Mad2l1bp', 'Maff', 'Malt1', 'Map1b', 'Map1lc3b', 'Mapk6', 'Marcks', 'Marcksl1', 'Mcl1', 'Med10', 'Mettl4', 'Mitf', 'Mllt11', 'Mmp13', 'Mndal', 'Mrpl52', 'Msantd3', 'Msmo1', 'Mt1', 'Mt2', 'Mthfd2', 'Mtmr7', 'Mvd', 'Nabp1', 'Narf', 'Nars', 'Nfe2l1', 'Nfil3', 'Nfkb1', 'Nfkb2', 'Nfkbia', 'Nfkbib', 'Nfkbid', 'Nfkbie', 'Nfkbiz', 'Ninj1', 'Nlrp3', 'Npnt', 'Nqo1', 'Nr4a3', 'Nsdhl', 'Nub1', 'Nupr1', 'Oasl1', 'Odc1', 'Oser1', 'Osgin1', 'Osgin2', 'P2rx4', 'Parp3', 'Pcdh15', 'Pdcd1', 'Pde4b', 'Pdgfb', 'Pdlim7', 'Pdp2', 'Pea15a', 'Pef1', 'Peli1', 'Pgd', 'Phlda1', 'Pik3r5', 'Pim1', 'Pinx1', 'Pla2g4e', 'Plagl2', 'Plaur', 'Plek', 'Plekhm3', 'Plekho2', 'Plin2', 'Plk2', 'Plk3', 'Plscr1', 'Pmaip1', 'Pmvk', 'Ppp1r11', 'Ppp1r15a', 'Ppp4r2', 'Prdx1', 'Prdx5', 'Prdx6', 'Prkg2', 'Procr', 'Prr13', 'Psmb6', 'Psmc2', 'Psmc6', 'Psmd10', 'Psmd11', 'Psmd7', 'Psmd8', 'Psph', 'Ptger2', 'Ptgr1', 'Ptgs2', 'Ptpn14', 'Ptpre', 'Pvr', 'RP23-172K2.12', 'RP23-57A17.2', 'Rab11fip1', 'Rab32', 'Rab8b', 'Rap1b', 'Rap2b', 'Rasgef1b', 'Rassf4', 'Rc3h1', 'Rdh11', 'Rel', 'Relb', 'Rffl', 'Rgs1', 'Rgs16', 'Rhbdf2', 'Rhob', 'Rhoc', 'Rhof', 'Riok3', 'Rnf19b', 'Rrs1', 'Saa3', 'Sars', 'Sat1', 'Sc5d', 'Sdc4', 'Sde2', 'Sele', 'Selk', 'Sept11', 'Serpinb1b', 'Serpinb2', 'Serpinb9b', 'Serpinb9e', 'Serpine1', 'Sertad1', 'Sertad2', 'Sgk1', 'Skil', 'Slc11a1', 'Slc15a3', 'Slc17a6', 'Slc1a4', 'Slc25a33', 'Slc25a37', 'Slc2a6', 'Slc39a1', 'Slc3a2', 'Slc48a1', 'Slc4a7', 'Slc6a9', 'Slc7a11', 'Slc7a6os', 'Slc9a4', 'Slfn10-ps', 'Slfn2', 'Slpi', 'Snx10', 'Socs3', 'Socs5', 'Sod2', 'Spp1', 'Spryd7', 'Spty2d1', 'Sqrdl', 'Sqstm1', 'Srfbp1', 'Srgn', 'Srxn1', 'Stap1', 'Stard4', 'Stip1', 'Stk40', 'Stx11', 'Stx6', 'Taf13', 'Taf7', 'Tank', 'Tars', 'Tbpl1', 'Tbrg1', 'Tfec', 'Tgif1', 'Tgoln1', 'Tiparp', 'Tlcd2', 'Tlr2', 'Tma16', 'Tmbim6', 'Tmem171', 'Tmem41b', 'Tmsb10', 'Tnf', 'Tnfaip2', 'Tnfaip3', 'Tnfaip8l2', 'Tnfrsf10b', 'Tnfrsf1b', 'Tnfsf4', 'Tnfsf9', 'Tnip1', 'Tnip3', 'Top1', 'Tox4', 'Tpbg', 'Tpm4', 'Traf1', 'Treml4', 'Trex1', 'Trib1', 'Trim13', 'Trpv2', 'Tslp', 'Tubb2a', 'Txndc9', 'Txnrd1', 'Uba6', 'Ubb', 'Ubc', 'Ube2h', 'Ubxn4', 'Unc80', 'Vegfa', 'Vmn2r90', 'Wsb1', 'Ypel5', 'Zbtb32', 'Zc3h12a', 'Zc3h12c', 'Zc3hav1', 'Zfand2a', 'Zfand5', 'Zfp330', 'Zfp36', 'Zfp617', 'Zfp622', 'Znrf1', 'Zswim4', 'Zwint', 'Zyx']
# list_of_genes = ['Cebpb', 'Isg15', 'Slc25a37', 'Ifih1', 'Peli1']

list_of_genes = ['Ccl2', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl6', 'Ccl7', 'Ccl9', 'Ccl11', 'Ccl17', 'Ccl20', 'Ccl22', 'Ccl24', 'Ccl25', 'Ccl28']
list_of_genes += ['Cxcl1', 'Cxcl2', 'Cxcl3', 'Cxcl5', 'Cxcl9', 'Cxcl10', 'Cxcl11', 'Cxcl14', 'Cxcl15', 'Cxcl16', 'Cxcl17']
# list_of_genes = list(set(list_of_genes) & set(list_of_all_genes))

print len(list_of_all_genes)
print list_of_genes

file_name = "trial_23_all_ccl_and_cxcl_genes.pdf"

"""
Analyze all the time points
"""

times_to_analyze = [0, 75, 150, 300]
cluster_list = {}
cluster_name_dict = {'0':{}, '75':{}, '150':{}, '300':{}}

for time_point in times_to_analyze:

	print "Analyzing " + str(time_point) + " minute time point"
	all_cells = []
	cell_names = []
	longest_time = 0
	number_of_cells = 0

	for cell in all_cells_total:
		if cell.time_point == time_point and cell.condition == 'Stim':
			number_of_cells += 1
			longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])
			all_cells += [cell]
			cell_names += [cell.id]

	dynamics_matrix = np.zeros((number_of_cells,longest_time), dtype = 'float32')

	"""
	Fill up the dynamics heat map matrix
	"""
	cell_counter = 0
	for cell in all_cells:
		dynam = cell.NFkB_dynamics
		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		cell_counter += 1

	"""
	Perform hierarchical clustering of the dynamics
	"""
	distance_matrix_dynamics = np.zeros((number_of_cells, number_of_cells))
	if time_point != 0:
		dynamics_load = np.load('/home/dvanva/SingleCellSequencing/' + str(time_point)+'_dynamics_distance_matrix_kshape.npz')
		distance_matrix_dynamics = dynamics_load["distance_matrix"]

		Y_dynamics = sch.linkage(distance_matrix_dynamics, method = 'ward')
		ind_dynamics = sch.fcluster(Y_dynamics,0.5*np.amax(Y_dynamics[:,2]),'distance')

	if time_point == 0:
		cluster_list[str(time_point)] = np.arange(1,2)
	else:
		cluster_list[str(time_point)] = np.arange(np.amin(ind_dynamics), np.amax(ind_dynamics)+1)


	if time_point == 0:
		for j in xrange(number_of_cells):
			all_cells[j].clusterID = 1
	else:
		for j in xrange(number_of_cells):
			all_cells[j].clusterID = ind_dynamics[j]

	cluster_dict = {}

	for cell in all_cells:
		cluster_dict[cell.id] = str(cell.clusterID)

	for cluster in cluster_list[str(time_point)]:
		cluster_name_dict[str(time_point)][str(cluster)] = []
		for cell in all_cells:
			if cell.clusterID == cluster:
				cluster_name_dict[str(time_point)][str(cluster)] += [cell.id]

"""
Compute posterior FPM distribution for a given gene
"""

plt.clf()
# fig, axes = plt.subplots(len(list_of_genes)/20+1,20, figsize = (4*20,4*(len(list_of_genes)/20 + 1)))
# fig, axes = plt.subplots(len(list_of_genes)/5+1,5, figsize = (4*5,4*(len(list_of_genes)/5 + 1)))
fig, axes = plt.subplots(len(list_of_genes),1, figsize = (4*1,4*(len(list_of_genes))))

counter = 0

for gene in list_of_genes:
	print gene
	cluster_1_mean = []
	cluster_2_mean = []
	cluster_3_mean = []

	cluster_1_low = []
	cluster_2_low = []
	cluster_3_low = []

	cluster_1_high = []
	cluster_2_high = []
	cluster_3_high = []

	for time_point in times_to_analyze:

		gene_name = """'""" + gene + """'"""

		r("o.prior = scde.expression.prior(models = o.ifm, counts = counts_data_int, length.out = 400, max.value = 10, show.plot = FALSE )")
		r("""gene_counts = counts_data_int[c(""" + gene_name + ""","mt-Atp8"),]""")

		ratio_list = []
		post_list = []


		for cluster in cluster_list[str(time_point)]:
			if time_point == 0:
				list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(time_point)][str(cluster)])
				r("list_of_cells = " + list_of_cells_r.r_repr())
				r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
				r("prior = o.prior")
				r("jp_0 = joint_posterior[" + gene_name + ",]")
				r("jp_0 = t(jp_0)")
				r("ratio = scde:::calculate.ratio.posterior(jp_0, jp_0, prior = o.prior, n.cores = 2)")

				ratios = 10 ** np.float32(pandas2ri.ri2py(r("colnames(ratio)")))
				post = np.float32(pandas2ri.ri2py(r("ratio")))

				ratio_list += [ratios]
				post_list += [post]

			else:
				list_of_cells_r = ro.vectors.StrVector(cluster_name_dict[str(time_point)][str(cluster)])
				r("list_of_cells = " + list_of_cells_r.r_repr())
				r("""joint_posterior = scde.posteriors(models = o.ifm[list_of_cells,], gene_counts, o.prior, n.cores = 4)""")
				r("jp = joint_posterior[" + gene_name + ",]")
				r("jp = t(jp)")
				r("ratio = scde:::calculate.ratio.posterior(jp, jp_0, o.prior, n.cores = 2)")

				ratios = 10 ** np.float32(pandas2ri.ri2py(r("colnames(ratio)")))
				post = np.float32(pandas2ri.ri2py(r("ratio")))

				ratio_list += [ratios]
				post_list += [post]

		# Give the clusters the proper order
		if time_point == 150:
			ratio_list_new = []
			post_list_new = []

			ratio_list_new += [ratio_list[2]]
			ratio_list_new += [ratio_list[0]]
			ratio_list_new += [ratio_list[1]]
			post_list_new += [post_list[2]]
			post_list_new += [post_list[0]]
			post_list_new += [post_list[1]]
			
			ratio_list = ratio_list_new
			post_list = post_list_new

		if time_point == 0:
			ratio = ratio_list[0]
			post = post_list[0]

			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]

			cluster_1_mean += ratios
			cluster_2_mean += ratios
			cluster_3_mean += ratios

			cluster_1_low += err_low
			cluster_2_low += err_low
			cluster_3_low += err_low

			cluster_1_high += err_high
			cluster_2_high += err_high
			cluster_3_high += err_high

		if time_point == 75:
			ratio= ratio_list[0]
			post = post_list[0]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_1_mean += ratios
			cluster_3_mean += ratios
			cluster_1_low += err_low
			cluster_3_low += err_low
			cluster_1_high += err_high
			cluster_3_high += err_high

			ratio = ratio_list[1]
			post = post_list[1]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_2_mean += ratios
			cluster_2_low += err_low
			cluster_2_high += err_high

		if time_point == 150:
			ratio = ratio_list[0]
			post = post_list[0]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_1_mean += ratios
			cluster_1_low += err_low
			cluster_1_high += err_high

			ratio = ratio_list[1]
			post = post_list[1]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_2_mean += ratios
			cluster_2_low += err_low
			cluster_2_high += err_high

			ratio = ratio_list[2]
			post = post_list[2]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_3_mean += ratios
			cluster_3_low += err_low
			cluster_3_high += err_high

		if time_point == 300:
			ratio = ratio_list[0]
			post = post_list[0]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_1_mean += ratios
			cluster_1_low += err_low
			cluster_1_high += err_high

			ratio = ratio_list[1]
			post = post_list[1]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_2_mean += ratios
			cluster_2_low += err_low
			cluster_2_high += err_high

			ratio = ratio_list[2]
			post = post_list[2]
			ratios = [ratio[np.argmax(post)]]
			cumsum = np.cumsum(post)
			err_low = [ratio[zero_crossing(cumsum, 0.16)]]
			err_high = [ratio[zero_crossing(cumsum, 0.84)]]
			cluster_3_mean += ratios
			cluster_3_low += err_low
			cluster_3_high += err_high

	"""
	Plot posteriors
	"""

	colors = ['g', 'r', 'b']

	cluster_1_low = np.array(cluster_1_low)
	cluster_2_low = np.array(cluster_2_low)
	cluster_3_low = np.array(cluster_3_low)

	cluster_1_high = np.array(cluster_1_high)
	cluster_2_high = np.array(cluster_2_high)
	cluster_3_high = np.array(cluster_3_high)

	cluster_1_mean = np.array(cluster_1_mean)
	cluster_2_mean = np.array(cluster_2_mean)
	cluster_3_mean = np.array(cluster_3_mean)

	max_val = np.amax(np.array([np.amax(cluster_1_high), np.amax(cluster_2_high), np.amax(cluster_3_high)]))
	axes.flatten()[counter].plot(times_to_analyze, cluster_1_mean, color = colors[0], linewidth = 1, label = 'Cluster 1')
	axes.flatten()[counter].plot(times_to_analyze, cluster_2_mean, color = colors[1], linewidth = 1, label = 'Cluster 2')
	axes.flatten()[counter].plot(times_to_analyze, cluster_3_mean, color = colors[2], linewidth = 1, label = 'Cluster 3')

	axes.flatten()[counter].fill_between(times_to_analyze, cluster_1_low, cluster_1_high, alpha = 0.1, color = colors[0])
	axes.flatten()[counter].fill_between(times_to_analyze, cluster_2_low, cluster_2_high, alpha = 0.1, color = colors[1])
	axes.flatten()[counter].fill_between(times_to_analyze, cluster_3_low, cluster_3_high, alpha = 0.1, color = colors[2])

	# axes.flatten()[counter].set_xlabel('Time (minutes)', fontsize = 12)
	# axes.flatten()[counter].set_ylabel('Fold change', fontsize = 12)
	axes.flatten()[counter].set_title(gene, fontsize = 16)
	axes.flatten()[counter].set_ylim([0,np.ceil(1.05*max_val)])
	axes.flatten()[counter].set_yticks([0,np.ceil(1.05*max_val)])
	axes.flatten()[counter].set_xlim([0, np.amax(times_to_analyze)])
	axes.flatten()[counter].set_xticks(times_to_analyze)



	# print [np.abs(cluster_3_low - cluster_3_mean), np.abs(cluster_3_high - cluster_3_mean)]
	# max_val = np.amax(np.array([np.amax(cluster_1_high), np.amax(cluster_2_high), np.amax(cluster_3_high)]))
	# axes.flatten()[counter].errorbar(times_to_analyze, cluster_1_mean, yerr = [np.abs(cluster_1_low - cluster_1_mean), np.abs(cluster_1_high - cluster_1_mean)], fmt = '-o', color = colors[0], ecolor = colors[0], linewidth = 2, label = 'Cluster 1')
	# axes.flatten()[counter].errorbar(times_to_analyze, cluster_2_mean, yerr = [np.abs(cluster_2_low - cluster_2_mean), np.abs(cluster_2_high - cluster_2_mean)], fmt = '-o', color = colors[1], ecolor = colors[1], linewidth = 2, label = 'Cluster 2')
	# axes.flatten()[counter].errorbar(times_to_analyze, cluster_3_mean, yerr = [np.abs(cluster_3_low - cluster_3_mean), np.abs(cluster_3_high - cluster_3_mean)], fmt = '-o', color = colors[2], ecolor = colors[2], linewidth = 2, label = 'Cluster 3')

	# axes.flatten()[counter].set_xlabel('Time (minutes)', fontsize = 16)
	# axes.flatten()[counter].set_ylabel('Fold change', fontsize = 16)
	# axes.flatten()[counter].set_title(gene, fontsize = 16)
	# axes.flatten()[counter].set_ylim([0,np.ceil(1.05*max_val)])
	# axes.flatten()[counter].set_yticks([0,np.ceil(1.05*max_val)])
	# axes.flatten()[counter].set_xlim([0, 1.05*np.amax(times_to_analyze)])
	# axes.flatten()[counter].set_xticks(times_to_analyze)
	counter += 1

	fig.tight_layout()
	plt.savefig("plots/" + file_name)



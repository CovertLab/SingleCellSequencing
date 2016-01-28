"""
Analysis!

Cluster the time traces and then plot a heatmap for the dynamics traces
"""

"""
Import python packages
"""

import collections
import itertools
import os
import numpy
import scipy
import pandas as pd
import numpy as np
import matplotlib
import cPickle as pickle
# matplotlib.use("Agg")
import matplotlib.pyplot as plt

tpm = np.arange(1, 1e6, .1)
k_list = [1e3, 1e4]
labels = ["1000 mRNAs", "10000 mRNAs"]
fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111)
counter = 0
for k in k_list:
	p_drop = (1-tpm/1e6)**k
	ax.plot(tpm, p_drop, label = labels[counter])
	ax.set_xscale('log')
	ax.set_xlabel('TPM', fontsize = 12)
	ax.set_ylabel('Dropout probability', fontsize = 12)
	ax.set_title('Dropout probability vs expression', y= 1.05, fontsize = 14)
	counter += 1

ax.legend()
fig.tight_layout()

plt.show()

plt.savefig("plots/dropout_est.pdf")


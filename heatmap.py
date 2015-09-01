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
from matplotlib import pyplot as plt
from matplotlib import gridspec

def cleanAxis(ax):
	ax.set_frame_on(False)
	for label in ax.axes.get_xticklabels():
		label.set_visible(False)
	for label in ax.axes.get_yticklabels():
		label.set_visible(False)
	for tick in ax.axes.get_xticklines():
		tick.set_visible(False)
	for tick in ax.axes.get_yticklines():
		tick.set_visible(False)
	for spine in ax.spines.values():
		spine.set_visible(False)

def heatmap(ax, data, cmap = 'Reds'):
	ax.imshow(data)
	plt.set_cmap(cmap)

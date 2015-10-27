"""
Test scripts to align data with STAR
"""

import seq_functions
direc = "/scratch/PI/mcovert/dvanva/sequencing/library1"
seq_functions.make_directories(direc, subdirecs_to_make = ["alignedSTAR"])
seq_functions.run_STAR(direc,"1-10_S14")
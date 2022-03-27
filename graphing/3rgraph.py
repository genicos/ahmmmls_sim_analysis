import math
import sys
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
import os

analysis = open("gss_3R.csv").readlines()

analysis = analysis[1:]
for i in range(len(analysis)):
    analysis[i] = analysis[i].split(",")


loci = []
lnl = []
for line in analysis:
    loci.append(float(line[0]))
    lnl.append(float(line[2]))

peaks_loc = []
peaks_lnl = []

cutoff = 15
dis = 6

for i in range(dis, len(lnl)-dis, 1):
    if lnl[i] > cutoff:
        peak = True
        for j in range(i-dis, i+dis, 1):
            if lnl[i] < lnl[j]:
                peak = False
        if(peak):
            peaks_loc.append(loci[i])
            peaks_lnl.append(lnl[i])
        


fig, ax = plt.subplots()
ax.scatter(peaks_loc, peaks_lnl, color="r", zorder = 2)

ax.plot(loci, lnl, zorder = 1)

ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)

ax.set_ylabel("likelihood ratio")
ax.set_xlabel("Chromosome 3R")
plt.show()
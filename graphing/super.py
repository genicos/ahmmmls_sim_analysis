import math
import sys
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
import os

fullloc = plticker.MultipleLocator(base=0.01)
halfloc = plticker.MultipleLocator(base=0.005)

analysis = open("gss_super_dom", "r").readlines()

print(analysis[0])
analysis = analysis[1:]
for i in range(len(analysis)):
    analysis[i] = analysis[i].split("\t")




sel_to_i = {"0.005":0, "0.01":1, "0.02":2, "0.05":3}
selections = [0.005, 0.01, 0.02, 0.05]

m_to_i = {"0.05":0, "0.1":1, "0.2":2, "0.5":3}
m = [0.05,0.1, 0.2, 0.5]

generations_to_i = {"100":0, "200":1, "500":2, "1000":3}
generations = [100, 200, 500, 1000]
scat_gen = []

for i in range(4):
    for j in range(20):
        scat_gen.append(generations[i] - 20 + 2*j)

# sel 
# m   
# gens

estimated_sel = [[None]*4,[None]*4,[None]*4,[None]*4]

for i in range(4):
    for j in range(4):
        estimated_sel[i][j] = [[],[],[],[]]

for line in analysis:
    estimated_sel[sel_to_i[line[3]]][m_to_i[line[2]]][generations_to_i[line[1]]].append(1 - float(line[5]))

fig, ax = plt.subplots(4, 4, constrained_layout=True, sharex=True,sharey='row')

for i in range(4):
    for j in range(4):
        x = scat_gen
        y = []
        for k in range(4):
            y.extend(estimated_sel[i][j][k])
        
        ax[i,j].tick_params(axis='x', which='major', labelsize=7)
        
        ax[i,j].yaxis.set_label_coords(-0.3, 0.5)

        if (i == 3):
            ax[i,j].set_xticks([100,200,500,1000])
        if(i == 0):
            ax[i,j].yaxis.set_major_locator(halfloc)
        if(i == 1):
            ax[i,j].set_yticks([0,0.01,0.02])
        if(i>1):
            ax[i,j].yaxis.set_major_locator(fullloc)
        if (i == 0):
            ax[i,j].set_title("m = "+str(m[j]))
        if (j == 0):
            ax[i,j].set_ylabel("sel = "+str(selections[i]))
        
        """
        plt.scatter(x, y)
        plt.axhline(y=1 - selections[i], color='r', linestyle='-')
        plt.legend(title="m = "+str(1 - m[j])+", sel = "+str(selections[i]))
        """
        ax[i,j].scatter(x, y, s = 5)
        ax[i,j].grid()
        ax[i,j].axhline(y=selections[i], color='r', linestyle='-')
        
        

plt.show()
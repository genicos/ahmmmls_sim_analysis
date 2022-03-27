import math
import sys
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
import os

analysis = open("new_dom_add_comp", "r").readlines()

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



good = [[None]*4,[None]*4,[None]*4,[None]*4]

for i in range(4):
    for j in range(4):
        good[i][j] = [0,0,0,0]

for line in analysis:
    good[sel_to_i[line[3]]][m_to_i[line[2]]][generations_to_i[line[1]]] = int(line[4])


bad = [[None]*4,[None]*4,[None]*4,[None]*4]

for i in range(4):
    for j in range(4):
        bad[i][j] = [20 - good[i][j][0],20 - good[i][j][1],20 - good[i][j][2],20 - good[i][j][3]]





xlabels = ["100","200","500","1000"]
x = np.arange(len(xlabels))
print(x)
width = 0.35


# sel 
# m   
# gens

fig, axe = plt.subplots(4, 4, constrained_layout=True, sharex=True,sharey='row')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False)

for i in range(4):
    for j in range(4):
        ax = axe[i,j]
        print(good[i][j])

        ax.bar(x, good[i][j], width, color='g')
        ax.bar(x, bad[i][j], width,bottom=good[i][j], color='r')

        ax.set_yticks([])
        ax.set_xticks(x)

        if (i == 3):
            ax.set_xticklabels(xlabels)
        if (i == 0):
            ax.set_title("m = " + str(m[j]))
        if (j == 0):
            ax.set_ylabel("sel = "+str(selections[i]))


plt.show()


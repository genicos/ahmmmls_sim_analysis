import math
import sys
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
import os


analysis = open("new_two_site_super", "r").readlines()

analysis = analysis[1:]
for i in range(len(analysis)):
    analysis[i] = analysis[i].split("\t")


sel_to_i = {"('0.01', '0.01')":0, "('0.01', '0.05')":1, "('0.05', '0.05')":2}
selections = ["(0.01, 0.01)", "(0.01, 0.05)", "(0.05, 0.05)"]

dist_to_i = {"0.005":0, "0.01":1, "0.02":2, "0.05":3}
dist = [0.005,0.01, 0.02, 0.05]



#sel
#dist


good = [[0,0,0,0],[0,0,0,0],[0,0,0,0]]
bad  = [[0,0,0,0],[0,0,0,0],[0,0,0,0]]

for line in analysis:
    if float(line[3]) > float(line[4]) + 2.996:
        good[sel_to_i[line[1]]][dist_to_i[line[2]]] += 1
    else:
        bad[sel_to_i[line[1]]][dist_to_i[line[2]]] += 1

print(good)
print(bad)


xlabels = ["0.005","0.01","0.02","0.05"]
x = np.arange(len(xlabels))
width = 0.35

fig, axe = plt.subplots(1, 3, constrained_layout=True, sharex=True,sharey='row')
for j in range(3):
    ax = axe[j]

    ax.bar(x, good[j], width, color='g')
    ax.bar(x, bad[j], width,bottom=good[j], color='r')

    ax.set_yticks([])
    ax.set_xticks(x)

    if(j == 1):
        ax.set_xlabel("Distance between sites in morgans")
        
    ax.set_xticklabels(xlabels)
        
    ax.set_title("sel = " + str(selections[j]))


plt.show()





ts1 = [[None]*4,[None]*4,[None]*4]
ts2 = [[None]*4,[None]*4,[None]*4]

for i in range(3):
    for j in range(4):
        ts1[i][j] = []
        ts2[i][j] = []

for line in analysis:
    ts1[sel_to_i[line[1]]][dist_to_i[line[2]]].append(float(line[5]))
    ts2[sel_to_i[line[1]]][dist_to_i[line[2]]].append(float(line[6]))



fig, ax = plt.subplots(3, 4, constrained_layout=True, sharex='col',sharey='row')

for i in range(3):
    for j in range(4):


        
        if (i == 0):
            ax[i,j].set_title("dist = "+str(dist[j]))
        if (j == 0):
            ax[i,j].set_ylabel("sel = "+str(selections[i]))

        
        ax[i,j].set_ylim([0.195,0.235])
        ax[i,j].set_yticks([0.2,0.21,0.22,0.23])
        if (j == 0):
            ax[i,j].set_xlim([0.194,0.2025])
            ax[i,j].set_xticks([0.195,0.2])
            ax[i,j].axvline(x=0.1975, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.2025, color='r', linestyle='-', zorder=2)
        if (j == 1):
            ax[i,j].axvline(x=0.195, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.205, color='r', linestyle='-', zorder=2)
        if (j == 2):
            ax[i,j].set_xlim([0.184,0.196])
            ax[i,j].set_xticks([0.185,0.19,0.195])
            ax[i,j].axvline(x=0.19, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.21, color='r', linestyle='-', zorder=2)
        if (j == 3):
            ax[i,j].set_xlim([0.1675,0.1925])
            ax[i,j].set_xticks([0.17,0.18,0.19])
            ax[i,j].axvline(x=0.175, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.225, color='r', linestyle='-', zorder=2)

        ax[i,j].scatter(ts1[i][j], ts2[i][j], s = 5,zorder=3)
        ax[i,j].grid()

plt.show()



selt1 = [[None]*4,[None]*4,[None]*4]
selt2 = [[None]*4,[None]*4,[None]*4]
selo1 = [[None]*4,[None]*4,[None]*4]

for i in range(3):
    for j in range(4):
        selt1[i][j] = []
        selt2[i][j] = []
        selo1[i][j] = []

for line in analysis:
    selt1[sel_to_i[line[1]]][dist_to_i[line[2]]].append(1 - float(line[8]))
    selt2[sel_to_i[line[1]]][dist_to_i[line[2]]].append(1 - float(line[9]))
    selo1[sel_to_i[line[1]]][dist_to_i[line[2]]].append(1 - float(line[10]))

print(selt1)
print(selt2)
print(selo1)

scat_gen = [[],[],[]]
for i in range(3):
    for j in range(20):
        scat_gen[i].append(i - 0.25 + 0.5/20*j)


fig, ax = plt.subplots(3, 4, constrained_layout=True, sharex='col',sharey='row')

for i in range(3):
    for j in range(4):


        ax[i,j].set_xticks([])

        if (i == 0):
            ax[i,j].set_title("dist = "+str(dist[j]))
        if (j == 0):
            ax[i,j].set_ylabel("sel = "+str(selections[i]))

        if(i == 0):
            ax[i,j].set_ylim([0,0.022])
            ax[i,j].axhline(y=0.01, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.01, color='r', linestyle='-', zorder=2)
        if(i == 1):
            ax[i,j].set_ylim([0,0.08])
            ax[i,j].axhline(y=0.01, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.05, color='r', linestyle='-', zorder=2)
        if(i == 2):
            ax[i,j].axhline(y=0.05, color='r', linestyle='-', zorder=2)
            ax[i,j].axhline(y=0.05, color='r', linestyle='-', zorder=2)
        
        ax[i,j].scatter(scat_gen[0], selt1[i][j], s = 5,zorder=3)
        ax[i,j].scatter(scat_gen[0], selt2[i][j], s = 5,zorder=3)
        ax[i,j].scatter(scat_gen[0], selo1[i][j], s = 5,zorder=3)
        ax[i,j].grid()

plt.show()
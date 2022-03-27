import math
import sys
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.ticker as plticker
import os

gen_miss = open("gen_miss_super", 'r').readlines()
m_miss = open("m_miss_super", 'r').readlines()


gen_miss = gen_miss[1:]
m_miss = m_miss[1:]


for i in range(len(gen_miss)):
    gen_miss[i] = gen_miss[i].split("\t")

for i in range(len(m_miss)):
    m_miss[i] = m_miss[i].split("\t")




real_gen_to_i = {"100":0, "1000":1}
real_gen = [100,1000]
real_m_to_i = {"0.1":0, "0.5":1}
real_m = [0.1,0.5]

m_miss_to_i = {"0.01":0,"0.2":1,"0.4":2,"0.7":3}
m_miss_array = [0.01, 0.2, 0.4, 0.7]
gen_miss_to_i = {"50":0,"200":1,"500":2,"800":3}
gen_miss_array = [50,200,500,800]



#realgen
#realm
#miss


m_miss_data = [[None]*2, [None]*2]

for i in range(2):
    for j in range(2):
        m_miss_data[i][j] = [[],[],[],[]]

gen_miss_data = [[None]*2, [None]*2]

for i in range(2):
    for j in range(2):
        gen_miss_data[i][j] = [[],[],[],[]]

for line in m_miss:
    m_miss_data[real_gen_to_i[line[1]]][real_m_to_i[line[2]]][m_miss_to_i[line[3]]].append((float(line[5]), float(line[4]) > float(line[6]) + 1 ))

for line in gen_miss:
    gen_miss_data[real_gen_to_i[line[1]]][real_m_to_i[line[2]]][gen_miss_to_i[line[3]]].append((float(line[5]), float(line[4]) > float(line[6]) + 1 ))



gen_good = [[None]*2,[None]*2]

for i in range(2):
    for j in range(2):
        gen_good[i][j] = [0,0,0,0]
        for k in range(4):
            for l in range(20):
                if(gen_miss_data[i][j][k][l][1]):
                    gen_good[i][j][k] += 1


gen_bad = [[None]*2,[None]*2]

for i in range(2):
    for j in range(2):
        gen_bad[i][j] = [0,0,0,0]
        for k in range(4):
            gen_bad[i][j][k] = 20 - gen_good[i][j][k]




m_good = [[None]*2,[None]*2]

for i in range(2):
    for j in range(2):
        m_good[i][j] = [0,0,0,0]
        for k in range(4):
            for l in range(20):
                if(m_miss_data[i][j][k][l][1]):
                    m_good[i][j][k] += 1


m_bad = [[None]*2,[None]*2]

for i in range(2):
    for j in range(2):
        m_bad[i][j] = [0,0,0,0]
        for k in range(4):
            m_bad[i][j][k] = 20 - m_good[i][j][k]

gen_xlabels = ["50","200","500","800"]
gen_x = np.arange(len(gen_xlabels))
width = 0.35




fig, axe = plt.subplots(2, 2, constrained_layout=True, sharex=True,sharey='row')
for i in range(2):
    for j in range(2):
        ax = axe[i,j]

        ax.bar(gen_x, gen_good[i][j], width, color='g')
        ax.bar(gen_x, gen_bad[i][j], width,bottom=gen_good[i][j], color='r')

        ax.set_yticks([])
        ax.set_xticks(gen_x)

        if (i == 1):
            ax.set_xticklabels(gen_xlabels)
        if (i == 0):
            ax.set_title("m = " + str(real_m[j]))
        if (j == 0):
            ax.set_ylabel("gens: "+str(real_gen[i]))


plt.show()




m_xlabels = ["0.01","0.2","0.4","0.7"]
m_x = np.arange(len(m_xlabels))
width = 0.35

fig, axe = plt.subplots(2, 2, constrained_layout=True, sharex=True,sharey='row')
for i in range(2):
    for j in range(2):
        ax = axe[i,j]

        ax.bar(m_x, m_good[i][j], width, color='g')
        ax.bar(m_x, m_bad[i][j], width,bottom=m_good[i][j], color='r')

        ax.set_yticks([])
        ax.set_xticks(m_x)

        if (i == 1):
            ax.set_xticklabels(m_xlabels)
        if (i == 0):
            ax.set_title("m = " + str(real_m[j]))
        if (j == 0):
            ax.set_ylabel("gens: "+str(real_gen[i]))

plt.show()
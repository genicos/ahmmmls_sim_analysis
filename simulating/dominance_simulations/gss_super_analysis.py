import sys
import math
from os import listdir
from os.path import isfile, join

output_dir = "outputs"

OUTfiles = [f for f in listdir(output_dir) if isfile(join(output_dir, f)) and f[0] == 'O']

selection_array = ['0.005', '0.01', '0.02', '0.05']
m_array = ['0.05', '0.1', '0.2', '0.5']
gen_array = ['100', '200', '500', '1000']

info = {}

for f in OUTfiles:
    params = f.split('_')
    params[0] = params[0][3:]
    ident = tuple(params[0:2])
    if not (ident in info):
        info[ident] = [0,[],[],[],[],int(params[2]),int(params[3]),int(params[4])]
    info[ident][0] = info[ident][0] + 1
    
    OUT = open(output_dir+"/"+f,'r').readlines()
    first = True
    for i in range(len(OUT)):
        line = OUT[i]
        if line[0] == 'G':
            next_line = OUT[i+1]
            
            lnl = float(line.split(' ')[-1])
            
            strength = 0
            if first:
                strength = float(next_line.split('\t')[2].split(' ')[0].split(',')[2])
                info[ident][1].append(lnl)
                info[ident][2].append(strength)
            else:
                strength = float(next_line.split('\t')[2].split(' ')[0].split(',')[0])
                info[ident][3].append(lnl)
                info[ident][4].append(strength)
            first = False
     
print("count","gens","m","sel","lnl0","sel0","lnl1","sel1", sep='\t')

for entry in info:
    sel = selection_array[info[entry][7]]
    m = m_array[info[entry][6]]
    gens = gen_array[info[entry][5]]
    count = info[entry][0]
    
    for i in range(count):
        print(count,gens,m,sel,info[entry][1][i],info[entry][2][i],info[entry][3][i],info[entry][4][i],sep='\t')

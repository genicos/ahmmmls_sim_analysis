import sys

from os import listdir
from os.path import isfile, join

output_dir = "outputs"

OUTfiles = [f for f in listdir(output_dir) if isfile(join(output_dir, f)) and f[0] == 'O']



sel = "0.02"

m_array = ['0.1', '0.5']
gen_array = ['100', '1000']

m_miss_array = ['0.01','0.2','0.4','0.7']
gen_miss_array = ['50','200','500','800']



info = {}

for f in OUTfiles:
    params = f.split('_')
    params[0] = params[0][3:]
    ident = tuple(params[0:2])
    if not (ident in info):
        info[ident] = [0,[],[],[],[],int(params[3]),int(params[2]),int(params[4])]
    info[ident][0] = info[ident][0] + 1
    
    OUT = open(output_dir+"/"+f,'r').readlines()
    
    for i in range(len(OUT)):
        line = OUT[i]
        if line[0] == 'G':
            next_line = OUT[i+1]
            
            lnl = float(line.split(' ')[-1])
            
            strength = float(next_line.split('\t')[2].split(' ')[0].split(',')[2])
            info[ident][1].append(lnl)
            info[ident][2].append(strength)
        if line[0] == 'l':
            next_line = OUT[i+1]
            
            lnl = float(line.split(' ')[-1])
            
            strength = float(next_line.split(',')[-1])
            info[ident][3].append(lnl)
            info[ident][4].append(strength)
            

# print(info)

print("count","gens","m","gen_miss","lnl0","sel0","lnl1","sel1", sep='\t')

for entry in info:
    gen_miss = gen_miss_array[info[entry][7]]
    m = m_array[info[entry][6]]
    gens = gen_array[info[entry][5]]
    count = info[entry][0]
    
    lnl0_mean = 0
    lnl0_var = 0
    sel0_mean = 0
    sel0_var = 0
    lnl1_mean = 0
    lnl1_var = 0
    sel1_mean = 0
    sel1_var = 0

    for i in range(count):
        lnl0_mean += info[entry][1][i]
        sel0_mean += info[entry][2][i]
        lnl1_mean += info[entry][3][i]
        sel1_mean += info[entry][4][i]
    

    lnl0_mean /= count
    sel0_mean /= count
    lnl1_mean /= count
    sel1_mean /= count

    for i in range(count):
        print(count,gens,m,gen_miss,info[entry][1][i],info[entry][2][i],info[entry][3][i], info[entry][4][i],sep='\t')

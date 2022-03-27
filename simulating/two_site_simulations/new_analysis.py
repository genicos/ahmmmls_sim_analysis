import sys

from os import listdir
from os.path import isfile, join

output_dir = "outputs"



OUTfiles = [f for f in listdir(output_dir) if isfile(join(output_dir, f)) and f[0] == 'O']


selections_array = [('0.01','0.01'),('0.01','0.05'),('0.05','0.05')]

dist_array = ['0.005','0.01','0.02','0.05']

info = {}

for f in OUTfiles:
    params = f.split('_')
    params[0] = params[0][3:]
    ident = tuple(params[0:2])
    if not (ident in info):
        info[ident] = [0,[],[],[],[],[],[],[],[],int(params[2]),int(params[3])]
    info[ident][0] = info[ident][0] + 1
    
    OUT = open(output_dir+"/"+f,'r').readlines()
     
    for i in range(len(OUT)):
        line = OUT[i]
        if line[0] == 'G':
            
            line1 = OUT[i+1]
            line2 = OUT[i+2]
            
            lnl = float(line.split(' ')[-1])
            
            site1_loc = float(line1.split('\t')[1])
            site2_loc = float(line2.split('\t')[1])
            
            
            site1_sel = float(line1.split(',')[-1])
            site2_sel = float(line2.split(',')[-1])

            if(site2_loc < site1_loc):
                temp = site1_loc
                site1_loc = site2_loc
                site2_loc = temp

                temp = site1_sel
                site1_sel = site2_sel
                site2_sel = temp

            info[ident][1].append(lnl)
            info[ident][3].append(site1_loc)
            info[ident][4].append(site2_loc)
            info[ident][6].append(site1_sel)
            info[ident][7].append(site2_sel)
        
        if line[0] == 'l':
            line1 = OUT[i+1]
            
            lnl = float(line.split(' ')[-1])
            
            site_loc = float(line1.split('\t')[1])
            
            site_sel = float(line1.split(',')[-1])

            info[ident][2].append(lnl)
            info[ident][5].append(site_loc)
            info[ident][8].append(site_sel)
            
     
print("count","sel","dist","two_site_lnl","one_site_lnl","two_site1","two_site2","one_site","two_site_sel1","two_site_sel2","one_site_sel", sep='\t')

for entry in info:
    sel = selections_array[info[entry][9]]
    dist = dist_array[info[entry][10]]
    
    count = info[entry][0]
    
    for i in range(count):
        print(count,sel,dist,info[entry][1][i],info[entry][2][i],info[entry][3][i],info[entry][4][i],info[entry][5][i],info[entry][6][i],info[entry][7][i],info[entry][8][i],sep='\t')

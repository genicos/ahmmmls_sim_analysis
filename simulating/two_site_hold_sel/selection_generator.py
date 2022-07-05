import sys
import subprocess
from math import floor
from random import random
from time import sleep

print(sys.argv)

batch = int(sys.argv[1])
piece = int(sys.argv[2])
gens  = int(sys.argv[3])
m     = int(sys.argv[4])
dist  = int(sys.argv[5])
index = int(sys.argv[6])



strength = (.01,.01)

distance = [.005, .01, .02, .05][dist]

generations = ["100","200","500","1000"][gens]

prop = [.05,.1,.2,.5][m]


selection = open("src/selections/selection"+str(batch)+"_"+str(piece)+"_"+str(index), "w")
selections = ""

selection.write("S\tA\t0\t")

a = 1
b = 1 - strength[0]/2
c = 1 - strength[0]

site = str(0.2 - distance/2)
a = str(a)
b = str(b)
c = str(c)
        
selection.write(site + "\t" + a + "\t" + b + "\t" + c + "\n")
selections += site + "\t" + a + "\t" + b + "\t" + c + "\n"


selection.write("S\tA\t0\t")

a = 1
b = 1 - strength[1]/2
c = 1 - strength[1]

site = str(0.2 + distance/2)
a = str(a)
b = str(b)
c = str(c)
        
selection.write(site + "\t" + a + "\t" + b + "\t" + c + "\n")
selections += site + "\t" + a + "\t" + b + "\t" + c + "\n"


selection.close()

selection_write = open("outputs/selection"+str(batch)+"_"+str(piece)+"_"+str(index), "w")
selection_write.write(selections)
selection_write.close()





# generate output file

output_string = generations+"\t0\t100\t0\tselam_outputs/selam_output"+str(batch)+"_"+str(piece)+"_"+str(index)+"\n"
output = open("src/outputs/output"+str(batch)+"_"+str(piece)+"_"+str(index), "w")

output.write(output_string)
output.close()


#generate demography file
demo = open("src/demographies/demography"+str(batch)+"_"+str(piece)+"_"+str(index), "w")
demo.write("pop1\tpop2\tsex\t0\t1\n")
demo.write("0\t0\tA\t10000\t10000\n")
demo.write("0\ta0\tA\t"+str(prop)+"\t0\n")
demo.write("0\ta1\tA\t"+str(1-prop)+"\t0\n")
demo.close()

import sys
import subprocess
from math import floor
from random import random
from time import sleep

print(sys.argv)

batch = int(sys.argv[1])
piece = int(sys.argv[2])
sel   = int(sys.argv[3])
dist  = int(sys.argv[4])
index = int(sys.argv[5])


#10
strength = [(.001, .001),(.001, .005), (.001, .01), (.001, .05), (.005, .005), (.005,.01), (.005,.05) , (.01,.01), (.01,0.05),(.05,.05)][sel]

#4
distance = [0.005,0.01,0.02,0.05][dist]


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

output_string = "500\t0\t100\t0\tselam_outputs/selam_output"+str(batch)+"_"+str(piece)+"_"+str(index)+"\n"
output = open("src/outputs/output"+str(batch)+"_"+str(piece)+"_"+str(index), "w")

output.write(output_string)
output.close()



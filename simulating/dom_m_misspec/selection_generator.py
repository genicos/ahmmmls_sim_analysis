import sys
import subprocess
from math import floor
from random import random
from time import sleep

print(sys.argv)

batch = int(sys.argv[1])
piece = int(sys.argv[2])
m     = int(sys.argv[3])
gens  = int(sys.argv[4])
mm    = int(sys.argv[5])
index = int(sys.argv[6])


prop = [0.1,0.5][m]

generations = ['100','1000'][gens]


selection = open("src/selections/selection"+str(batch)+"_"+str(piece)+"_"+str(index), "w")
selections = ""


selection.write("S\tA\t0\t")


selection.write("0.2\t1\t1\t0.98\n")
selections += "0.2\t1\t1\t0.98\n"


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
demo = open("src/demographies/demography"+str(batch)+"_"+str(piece)+"_"+str(index),"w")
demo.write("pop1\tpop2\tsex\t0\t1\n")
demo.write("0\t0\tA\t10000\t10000\n")
demo.write("0\ta0\tA\t"+str(prop)+"\t0\n")
demo.write("0\ta1\tA\t"+str(1-prop)+"\t0\n")
demo.close()

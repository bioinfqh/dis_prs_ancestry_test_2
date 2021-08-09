#!/bin/bash

import os
import sys
cutoff = 0.01

percentages = {}
ctr = int(sys.argv[2])
for i in range(0,int(sys.argv[2])):
    file = open(str(sys.argv[1])+"."+str(i) +".Q")
    if(os.stat(str(sys.argv[1])+"."+str(i) +".Q").st_size == 0):
        ctr = ctr - 1
    for line in file:
        nat = line.split()[0]
        perc = float(line.split()[1])
        if(perc > cutoff):
            if(nat in percentages):
                percentages[nat] += float(perc)
            else:
                percentages[nat] = float(perc)

file.close()
outstr = ""
for nat in percentages:
    average = percentages[nat] / float(ctr)
    if(average > cutoff):
        outstr = outstr + nat + " " + str(average) + "\n"


with open("output.Q",'w') as myfile:
    myfile.write(outstr.rstrip())

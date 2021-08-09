#!/bin/bash

import sys

nat_dict = {}
conv_file = open(sys.argv[2])
for line in conv_file:
    nat_original = line.split()[0]
    nat_new = line.split()[1]
    nat_dict[nat_original] = nat_new

conv_file.close()

outstr = ""
pred_file = open(sys.argv[1])
for line in pred_file:
    id = line.split()[0]
    real_nat = line.split()[1]
    pred_nat = line.split()[2]
    if(real_nat in nat_dict):
        real_nat = nat_dict[real_nat]
    if(pred_nat in nat_dict):
        pred_nat = nat_dict[pred_nat]
    outstr = outstr + str(id) + "\t" + real_nat + "\t" + pred_nat + "\n"
pred_file.close()

print(outstr.rstrip())

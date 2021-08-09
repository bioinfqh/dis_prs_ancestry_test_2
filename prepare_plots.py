#!/bin/bash
import os
import numpy as np 
import pandas as pd
import sys

predictions = []
real_nationalities = []
list_of_nationalities = []
list_of_dicts = []

file = open(sys.argv[1])
ids = []
for line in file:
    id = line.split()[0]
    ids.append(id)
file.close()
for i in ids:
    file = open("validation/output_" + str(i) + ".txt")
    nats = []
    preds = []
    pred_dict = {}
    for line in file:
        nats.append(line.split()[0])
        preds.append(line.split()[1])
        #print(line.split()[0])
        if not (line.split()[0] in list_of_nationalities):
            list_of_nationalities.append(line.split()[0])
        pred_dict[line.split()[0]] = line.split()[1]
    nats_arr = np.asarray(preds,dtype=float)
    curr_max=np.argmax(nats_arr)
    #predictions.append(nats[curr_max])
    list_of_dicts.append(pred_dict)
    #outstr = outstr + "\n" + sample_id + "\t" + nationality + "\t" + nats[curr_max]
    #print("\n" + sample_id + "\t" + nationality + "\t" + nats[curr_max])
    file.close()
list_of_real_nationalities = []
file = open("validation/predicted_nationalities.txt")
for line in file:
    nationality = line.split()[1]
    real_nationalities.append(nationality)
    predictions.append(line.split()[2])
    if not (nationality in list_of_real_nationalities):
        list_of_real_nationalities.append(nationality)
#print(len(real_nationalities))
#print(len(predictions))
#print(list_of_nationalities)
predict_perc = pd.DataFrame(np.zeros((len(list_of_real_nationalities), len(list_of_nationalities))))
admix_perc = pd.DataFrame(np.zeros((len(list_of_real_nationalities), len(list_of_nationalities))))
predict_perc.index = list_of_real_nationalities
predict_perc.columns = list_of_nationalities
admix_perc.index = list_of_real_nationalities
admix_perc.columns = list_of_nationalities
#print(predict_perc)
#print(admix_perc)
for i in range(1,len(predictions)):
    predict_perc.at[real_nationalities[i],predictions[i]] += 1
    curr_dict = list_of_dicts[i]
    for key in curr_dict:
        admix_perc.at[real_nationalities[i],key] += float(curr_dict[key])

#print(predict_perc)
#print(admix_perc)
predict_perc.to_csv(r"validation/predicted_percentages.txt")
admix_perc.to_csv(r"validation/admixture_proportions.txt")
#for i in range(0:len(predictions)):
    


import sys
import operator


results = open(sys.argv[1])
patient_id = str(sys.argv[2])
outfile = str(sys.argv[3])
print_risks = "true"
id_col = 2
risk_col = 6
htmlstr=results.read()
results.close()
lines= htmlstr.split("\n")
scores = {}
for line in lines[1:]:
    lineSplit = line.split()
    #print(lineSplit)
    if(len(lineSplit) < 2):
        continue
    scores[str(lineSplit[id_col-1])] = float(lineSplit[risk_col-1])
    
sorted_preds_tmp = dict(sorted(scores.items(), key=operator.itemgetter(1)))
sorted_preds = [float(sorted_preds_tmp[i]) for i in sorted_preds_tmp]
sorted_ids = [key for key in sorted_preds_tmp]
#print(sorted_ids)
#print(sorted_preds[1])
#print(sorted_preds[100])
idx_of_patient = sorted_ids.index(patient_id)
relpos_pat = float(idx_of_patient) / float(len(sorted_ids))
if((float(sorted_preds[len(sorted_preds)-1]) - float(sorted_preds[0])) == 0):
    score_pat = 50.0
score_pat = (float(sorted_preds[idx_of_patient]) - float(sorted_preds[0])) / (float(sorted_preds[len(sorted_preds)-1]) - float(sorted_preds[0]))
#print(str(float(sorted_preds[idx_of_patient])))
#print(str(float(sorted_preds[0])))
#print(str(float(sorted_preds[len(sorted_preds)-1])))
if(print_risks == "true"):
    fh=open(outfile,'w')
    score_str_array = [(str(i) + "\t" + str(scores[i]) + "\t-9") for i in scores]
    fh.write("\n".join(score_str_array))
    fh.close()
print(patient_id + "\t" + str(score_pat) + "\t" + str(relpos_pat))

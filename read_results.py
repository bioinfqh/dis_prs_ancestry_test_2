
import sys


results = open(sys.argv[1])
target_id = str(sys.argv[2])
if(len(sys.argv) > 3):
    prs_pos = int(str(sys.argv[4])) - 1
    id_pos = int(str(sys.argv[3])) - 1
htmlstr=results.read()
lines= htmlstr.split("\n")
firstline=lines[0].split()
#print(firstline)
if(len(sys.argv) < 4):
    prs_pos = firstline.index("PRS")
    id_pos = firstline.index("IID")
prs_dict={}
for i in range(1,len(lines)):
    line = lines[i]
    if(len(line.split()) < prs_pos):
        continue
    name = line.split()[id_pos]
    prs = line.split()[prs_pos]
    prs_dict[str(name)]=float(prs)
#print(target_id)
#print(prs_dict)

ret_arr = [(key + ":" + str(prs_dict[key])) for key in prs_dict]
#print("\n".join(ret_arr))
values = [prs_dict[key] for key in prs_dict]
if(target_id in prs_dict):
    value_of_target = prs_dict[target_id]
    values_sorted = sorted(values,key=float)
    index = values_sorted.index(value_of_target)
    perc = (float(index) + 0.5) / float(len(values_sorted))
    print(str(perc))



#print(str(prs_dict[target_id]))

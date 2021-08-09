import sys
import math

file = open(sys.argv[1])
col = int(str(sys.argv[2])) - 1
#target_id = str(sys.argv[2])

htmlstr=file.read()
lines= htmlstr.split("\n")
linesRet = []
lines = lines[1:]
for line in lines:
    lineSplit = line.split()
    if(len(lineSplit) < col):
        continue
    pval = float(lineSplit[col])
    lineSplit[col] = str(math.log(pval))
    linesRet.append("\t".join(lineSplit))

print("\n".join(linesRet))

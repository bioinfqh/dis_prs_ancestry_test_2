import sys

path = sys.argv[1]
table_file = open(path)
htmlstr=table_file.read()
table_file.close()
lines = htmlstr.split("\n")
for line in lines:
    if("#CHROM") in line:
        lineSplit = line.split()
        sample_id = str(lineSplit[len(lineSplit)-1])
        print(sample_id)
        sys.exit()
    

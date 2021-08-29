import json



def read_anc(path):
    table_file = open(path)
    htmlstr=table_file.read()
    table_file.close()
    lines = htmlstr.split("\n")
    anc_dict = {}
    for line in lines:
        if(len(line.split()) < 2):
            continue
        nat = line.split()[0]
        perc = line.split()[1]
        if("." in perc):
            anc_dict[nat] = float(perc)
    return(anc_dict)


def anc_to_json(anc_dict):
    ret_str = ""
    json_str = json.dumps(anc_dict)
    return(json_str)


result_file = str(sys.argv[1])
outfile = str(sys.argv[2])

anc_dict = read_anc(result_file)
#dict_new = {}
countries_and_percentages = []
for nat in anc_dict:
    dict_new = {}
    dict_new['nat'] = nat
    dict_new['perc'] = anc_dict[nat]
    countries_and_percentages.append(dict_new)
anc_json = anc_to_json(anc_dict)
fh=open(outfile,'w')
fh.write(anc_json)
fh.close()

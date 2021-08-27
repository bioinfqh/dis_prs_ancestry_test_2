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

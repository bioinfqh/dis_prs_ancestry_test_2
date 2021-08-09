import os
from dis_calc.extract_table import run_all
#from extract_table import run_all
import pandas as pd
import json

get_abs_risk = "f"
disease_group_file = "dis_calc/disease_groups.txt"


## this method parses the sample name from a vcf file (last tab of the header line)
def get_id_from_vcf(path):
    table_file = open(path)
    htmlstr=table_file.read()
    table_file.close()
    lines = htmlstr.split("\n")
    for line in lines:
        if("#CHROM") in line:
            lineSplit = line.split()
            sample_id = str(lineSplit[len(lineSplit)-1])
            return(sample_id)
    return("none")

## this method reads the output.Q file (result from ancestry) into a dict
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
        
## PRS calculation. patient_id and customer_id are usually identical
def make_prs(bg_dataset,disease_list_path,patient_vcf,patient_id,customer_id):
    # get patient ID
    patient_id_new = get_id_from_vcf(str(patient_vcf))
    # run PRS wrapper script
    bash_code = 'bash dis_calc/prs_and_report.sh '+ str(bg_dataset) + ' ' + str(patient_vcf) + ' ' + str(disease_list_path) + ' ' + str(patient_id_new) + ' ' + str(get_abs_risk)
    #bash_code = 'bash dis_calc/prs_and_report.sh '+ str(bg_dataset) + ' ' + str(patient_vcf) + ' ' + str(disease_list_path) + ' ' + str(patient_id) + ' ' + str(get_abs_risk)
    os.system(bash_code)
    # copy result to static directory
    report_file_name = "dis_calc/static/prs_report_" + customer_id + ".pdf"
    prs_json_path = "dis_calc/static/prs_" + patient_id_new + ".json"
    os.system('mv your_prs_report.pdf ' + report_file_name)
    os.system('mv ' + prs_json_path + ' dis_calc/static/prs_' + customer_id + '.json')
    
## disease gene calculation. patient_id and customer_id are usually the same
def make_dis_gene(patient_vcf,patient_id,dis_group_file,customer_id):
    diseases = []
    diseases_to_ret = []
    if(dis_group_file == "all"):
        diseases = ["all"]
        diseases_to_ret = ["all"]
    else:
        # read list of diseases from file
        disgr_file = open(dis_group_file)
        disgr_str=disgr_file.read()
        disgr_file.close()
        lines = disgr_str.split("\n")
        for line in lines:
            lineSplit = line.split("\t")
            disgr_name = lineSplit[0].replace(" ","_")
            diseases.append(disgr_name)
            diseases_to_ret.append(lineSplit[0])
            #diseases = diseases + "," + disgr_name
    vcf_df = pd.read_csv(patient_vcf,sep='\t')
    # run analysis (from make_prs_report)
    report_paths = run_all(vcf_df,diseases,("dis_calc/static/dis_report_" + customer_id))
    # remove multiple static-paths in link
    for i in report_paths:
        i = i.replace("dis_calc/static/","")
    return(report_paths,diseases_to_ret)
    #bash_code = 'python3 extract_table.py --input ' + patient_vcf + ' --disease ' + diseases + ' --output static/your_report_' + customer_id
    #os.system(bash_code)

## ancestry estimation.
def make_ancestry_report(patient_vcf,patient_id):
    outfile = "dis_calc/static/ancestry_map_" + patient_id + ".html"
    outfile_for_ret = "ancestry_map_" + patient_id + ".html"
    result_file = "output_" + patient_id + ".Q"
    # run ancestry wrapper script
    bash_code = 'bash get_ancestry.sh dis_calc/anc_vcf.vcf ' + outfile + ' ' + result_file
    os.system(bash_code)
    anc_dict = read_anc(result_file)
    #dict_new = {}
    countries_and_percentages = []
    for nat in anc_dict:
        dict_new = {}
        dict_new['nat'] = nat
        dict_new['perc'] = anc_dict[nat]
        countries_and_percentages.append(dict_new)
    anc_json = anc_to_json(anc_dict)
    fh=open("dis_calc/static/ancestry_" + patient_id + ".json",'w')
    fh.write(anc_json)
    fh.close()
    return(countries_and_percentages,outfile_for_ret,result_file)

def make_local_anc_report(patient_vcf,customer_id):
    anc_img = "chm_img_" + customer_id + ".png"
    sample_id =get_id_from_vcf(str(patient_vcf))
    bash_code = 'bash run_local_anc.sh ' + patient_vcf + ' ' + sample_id + ' ' + customer_id
    os.system(bash_code)
    return(anc_img)


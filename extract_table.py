
import sys
import urllib.request
import urllib.error
import urllib.parse
import pandas as pd
import numpy as np
from rsid_to_clinvar import get_hgvs_list
from rsid_to_clinvar import get_correct_hgvs_link
from rsid_to_clinvar import get_correct_hgvs
from rsid_to_clinvar import get_inheritance_mode
from rsid_to_clinvar import get_publications
from rsid_to_clinvar import get_synonyms
from rsid_to_clinvar import get_ticket
from rsid_to_clinvar import is_synonym
from rsid_to_clinvar import get_inheritance_mode_2
from rsid_to_clinvar import get_inheritance_mode_3
from rsid_to_clinvar import get_link
from rsid_to_clinvar import read_gene_associations
from rsid_to_clinvar import read_gene_associations_2
from rsid_to_clinvar import read_pop_dict
from rsid_to_clinvar import get_synonyms_new
#import pdfkit 
import mygene
import json
from requests.auth import HTTPBasicAuth
import requests
import os

  

run_as_script = "false"

# to always prefer pathogenic variants: true, if choice of laboratory more important: false
sort_by_clin_sig = "false"


path_prefix = "/scripts/"

assoc_table_path = str(path_prefix + "tableExport.csv")
gene_disease_groups_path = str(path_prefix + "gene_disease_groups.csv")

url = 'https://clinvarminer.genetics.utah.edu/submissions-by-variant/NM_001354689.3%28RAF1%29%3Ac.1142G%3EC%20%28p.Gly381Ala%29'
url = 'https://clinvarminer.genetics.utah.edu/submissions-by-variant/NM_000059.4(BRCA2)%3Ac.4258G%3ET%20(p.Asp1420Tyr)'

mg = mygene.MyGeneInfo()

ret_all = {}
gene_name_dict = {}
prot_change_dict = {}
dict_by_snpid = {}
snp_data = {}

polyphen_threshold = 0.9
remove_benign = "false"


synlist_path = str(path_prefix + "syndict_temp.txt")

run_as_script="true"
use_all_lines="false"
automatic_stop_line_def="true"

start_line_def = 6600
stop_line_def = 6700

pubmed_df_max_nr = 4

pop_dict={}
pop_dict["AFR"] = "African"
pop_dict["AMR"] = "Mixed American"
pop_dict["ASJ"] = "Ashkenazi Jewish"
pop_dict["EAS"] = "East Asian"
pop_dict["FIN"] = "Finnish in Finland"
pop_dict["NFE"] = "Non-Finnish European"
pop_dict["SAS"] = "South Asian"
pop_dict=read_pop_dict(str(path_prefix + "pop_list.txt"))


file_with_labranks = str(path_prefix + "list_of_labranks.txt")
file_with_varsome_terms = str(path_prefix +"list_of_terms_for_varsome.txt")
#def make_pdf(html_path,outfile):
#    pdfkit.from_file(html_path, outfile,options={'page-height': '317mm', 'page-width': '210mm' , 'dpi':400})


def write_json_to_files(strings,outfiles):
    if not(len(strings) == len(outfiles)):
        return("none")
    for i in range(0,len(strings)):
        curr_string = strings[i]
        curr_outfile = outfiles[i]
        fh=open(curr_outfile,'w')
        fh.write(curr_string)
        fh.close()
    return(outfiles)

def write_json_to_one_file(strings,outfile):
    #print(strings)
    fh=open(outfile,'w')
    for i in range(0,len(strings)):
        curr_string = strings[i]
        fh.write(curr_string)
        if not((i == 0) and (len(strings) == 1)):
            fh.write("\n")
    fh.close()
    return(outfile)
        
def generate_json_OLD(result_text,bgcolor,ret_df,comment_df,pubmed_df):
    df_for_html = ret_df[["clin_sig","gene_name","exon","rsid","dna_change","prot_change","transcript","zyg","inh","diseases_STR"]]
    df_for_html.columns=["Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Disease(s)"]
    drug_response_df = df_for_html[df_for_html['Classification']=="drug response"]
    new_col = np.arange(1,(len(df_for_html) + 1))
    col_nbr_dict =  dict(zip(df_for_html["SNP ID"].tolist(),new_col))
    df_for_html.insert(0,"Number",new_col)
    drug_response_df.columns=["Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Response(s)"]
    df_for_html_2 = df_for_html[df_for_html['Classification']!="drug response"]
    df_for_html = df_for_html_2
    #print(df_for_html)
    #print(df_for_html["Classification"].tolist())
    #df_for_html['nbr'] = np.arange(len(df_for_html))
    #df_for_html = df_for_html[["nbr","clin_sig","gene_name","exon","rsid","dna_change","prot_change","transcript","zyg","inh","diseases_STR"]]
    #df_for_html.columns=["Number","Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Disease(s)"]
    #print(list(ret_df.columns))
    df_for_freq = ret_df[["clin_sig","gene_name","transcript","dna_change","prot_change","chr_loc","freq","rsid"]]
    df_for_freq = df_for_freq[df_for_freq['clin_sig']!="drug response"]
    nbr_col_2 = [col_nbr_dict[rsid] for rsid in df_for_freq["rsid"].tolist()]
    df_for_freq = df_for_freq[["gene_name","transcript","dna_change","prot_change","chr_loc","freq","rsid"]]
    df_for_freq.columns=["Gene","Transcript","DNA change","Protein change","Genomic Location","Alternate Allele Fraction","dbSNP rsID"]
    #print(df_for_freq)
    #print(df_for_freq["rsid"].tolist())
    #print(nbr_col_2)
    df_for_freq.insert(0,"Mutation\nNumber",nbr_col_2)
    comment_df = comment_df[comment_df['clin_sig']!="drug response"]
    comment_df = comment_df[comment_df['rsid'].notna()]
    nbr_col_3 = [col_nbr_dict[rsid] for rsid in comment_df["rsid"].tolist()]
    #print(comment_df["clin_sig"].tolist())
    comment_df.insert(0,"Number",nbr_col_3)
    if not(len(list(pubmed_df.index)) < 1):
        nbr_col_4 = [col_nbr_dict[rsid] for rsid in pubmed_df['rsid'].tolist()]
        pubmed_df.insert(0,"Number",nbr_col_4)
        pubmed_df = pubmed_df[["Number","name","PMID","title"]]
        pubmed_df.columns = ["Number","Name of Mutation","Pubmed ID","Title"]
    results = df_for_html.to_dict(orient='records')
    comments = comment_df.to_dict(orient='records')
    freqs = df_for_freq.to_dict(orient='records')
    pubmed = pubmed_df.to_dict(orient='records')
    results_json = json.dumps(results)
    comments_json = json.dumps(comments)
    freqs_json = json.dumps(freqs)
    pubmed_json = json.dumps(pubmed)
    return(results_json,comments_json,freqs_json,pubmed_json)

def generate_json(result_text,bgcolor,ret_df,comment_df,pubmed_df,cancer_groups_filter_active,cancer_groups_for_filter):
    df_for_html = ret_df[["clin_sig","gene_name","exon","rsid","dna_change","prot_change","transcript","zyg","inh","diseases_STR","cancer_groups"]]
    df_for_html.columns=["Classification","Gene","Exon/Intron","SNP ID","DNA change","Protein change","Transcript ID","Zygosity","Inheritance","Associated Disease(s)","Cancer Groups"]
    drug_response_df = df_for_html[df_for_html['Classification']=="drug response"]
    new_col = np.arange(1,(len(df_for_html) + 1))
    col_nbr_dict =  dict(zip(df_for_html["SNP ID"].tolist(),new_col))
    df_for_html.insert(0,"Number",new_col)
    drug_response_df.columns=["Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Response(s)","Cancer Groups"]
    df_for_html_2 = df_for_html[df_for_html['Classification']!="drug response"]
    df_for_html = df_for_html_2
    #print(df_for_html)
    #print(df_for_html["Classification"].tolist())
    #df_for_html['nbr'] = np.arange(len(df_for_html))
    #df_for_html = df_for_html[["nbr","clin_sig","gene_name","exon","rsid","dna_change","prot_change","transcript","zyg","inh","diseases_STR"]]
    #df_for_html.columns=["Number","Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Disease(s)"]
    #print(list(ret_df.columns))
    df_for_freq = ret_df[["clin_sig","gene_name","transcript","dna_change","prot_change","chr_loc","freq","rsid"]]
    df_for_freq = df_for_freq[df_for_freq['clin_sig']!="drug response"]
    nbr_col_2 = [col_nbr_dict[rsid] for rsid in df_for_freq["rsid"].tolist()]
    df_for_freq = df_for_freq[["gene_name","transcript","dna_change","prot_change","chr_loc","freq","rsid"]]
    df_for_freq.columns=["Gene","Transcript","DNA change","Protein change","Genomic Location","Alternate Allele Fraction","dbSNP rsID"]
    #print(df_for_freq)
    #print(df_for_freq["rsid"].tolist())
    #print(nbr_col_2)
    df_for_freq.insert(0,"Mutation\nNumber",nbr_col_2)
    comment_df = comment_df[comment_df['clin_sig']!="drug response"]
    comment_df = comment_df[comment_df['rsid'].notna()]
    nbr_col_3 = [col_nbr_dict[rsid] for rsid in comment_df["rsid"].tolist()]
    #print(comment_df["clin_sig"].tolist())
    comment_df.insert(0,"Number",nbr_col_3)
    if not(len(list(pubmed_df.index)) < 1):
        nbr_col_4 = [col_nbr_dict[rsid] for rsid in pubmed_df['rsid'].tolist()]
        pubmed_df.insert(0,"Number",nbr_col_4)
        pubmed_df = pubmed_df[["Number","name","PMID","title"]]
        pubmed_df.columns = ["Number","Name of Mutation","Pubmed ID","Title"]
    comments = comment_df.to_dict(orient='records')
    freqs = df_for_freq.to_dict(orient='records')
    pubmed = pubmed_df.to_dict(orient='records')
    #print(comments)
    #print(freqs)
    #print(pubmed)
    result_dict = df_for_html.to_dict(orient="records")
    dict_1 = {}
    ctr = 0
    cancer_groups_filter_active = "true"
    cancer_groups_for_filter = ["LUNG","ENDOCRINE"]
    dict_all = {}
    for gene_name_curr in list(df_for_html.index):
        curr_nbr_gene = str(df_for_html.loc[gene_name_curr,"Number"])
        gene_name_tmp = str(df_for_html.loc[gene_name_curr,"Gene"])
        gene_loc_tmp = str(df_for_html.loc[gene_name_curr,"DNA change"])
        gene_name = gene_name_tmp + " " + gene_loc_tmp
        new_dict = {}
        list_of_comments = []
        list_of_freqs = []
        list_of_pubmeds = []
        curr_row_as_dict = {}
        genes_found_by_cancer_type = {}
        for i in range(0,len(cancer_groups_for_filter)):
            genes_found_by_cancer_type[str(cancer_groups_for_filter[i])] = []
        if(cancer_groups_filter_active == "true"):
            curr_gene_cancer_groups = df_for_html.loc[gene_name_curr,"Cancer Groups"]
            for tmp in cancer_groups_for_filter:
                if(tmp in curr_gene_cancer_groups):
                    genes_found_by_cancer_type[tmp].append(gene_name)
            #if(found == "false"):
            #    continue
        if not(len(result_dict) < ctr):
            result_row = result_dict[ctr]
            ctr = ctr + 1
            diseases_for_curr_row = str(result_row["Associated Disease(s)"])
            disease_list_tmp = diseases_for_curr_row.split(";")
            s = []
            s_new = []
            for i in disease_list_tmp:
                #print(i)
                i_new = str("".join(i.split())).replace(" ","").replace("\t", "").lstrip()
                #print(i_new)
                #i_new = i_new.replace(" ","").replace("\t", "")
                if i_new.lower() not in s_new:
                    s.append(i)
                    s_new.append(i_new.lower())
                #print(s)
                #print(s_new)
            diseases_current = s
            disease_list_tmp = diseases_current
            result_row["Associated Disease(s)"] = disease_list_tmp
            curr_row_as_dict = result_row
        for i in range(0,len(comments)):
            if not('Number' in comments[i]):
                continue
            curr_nbr = str(comments[i]['Number'])
            if(curr_nbr == curr_nbr_gene):
                del comments[i]["Number"]
                del comments[i]["gene_name"]
                del comments[i]["rsid"]
                tmp_dict = comments[i]
                list_of_comments.append(tmp_dict)
        for i in range(0,len(freqs)):
            if not('Mutation\nNumber' in freqs[i]):
                continue
            curr_nbr = str(freqs[i]['Mutation\nNumber'])
            if(curr_nbr == curr_nbr_gene):
                del freqs[i]["Mutation\nNumber"]
                tmp_dict = freqs[i]
                list_of_freqs.append(tmp_dict)
                frec_as_dict = tmp_dict
        for i in range(0,len(pubmed)):
            if not('Number' in pubmed[i]):
                continue
            curr_nbr = str(pubmed[i]['Number'])
            if(curr_nbr == curr_nbr_gene):
                del pubmed[i]["Number"]
                del pubmed[i]["Name of Mutation"]
                tmp_dict = pubmed[i]
                list_of_pubmeds.append(tmp_dict)
        new_dict['Gene Data'] = curr_row_as_dict
        new_dict['Comments'] = list_of_comments
        new_dict['Frequencies'] = frec_as_dict
        new_dict['Publications'] = list_of_pubmeds
        dict_1[gene_name] = new_dict
        if(cancer_groups_filter_active == "true"):
            tmp = {}
            for group in cancer_groups_for_filter:
                if(gene_name in genes_found_by_cancer_type[group]):
                    tmp[gene_name] = new_dict
                    if not(group in dict_all):
                        dict_all[group] = [tmp]
                    else:
                        dict_all[group].append(tmp)
    #results = df_for_html.to_dict(orient='records')
    if(cancer_groups_filter_active == "true"):
        results_json = json.dumps(dict_all)
    else:
        results_json = json.dumps(dict_1)
    results_json = results_json.replace("\n"," ").replace("NaN","\"NaN\"").replace("\"\"NaN\"\"","\"NaN\"")
    comments_json = json.dumps(comments)
    freqs_json = json.dumps(freqs)
    pubmed_json = json.dumps(pubmed)
    return(results_json,comments_json,freqs_json,pubmed_json)


def generate_json_2(result_text,bgcolor,ret_df,comment_df,pubmed_df,cancer_groups_filter_active,cancer_groups_for_filter,customer_id,varsome_dict):
    df_for_html = ret_df[["clin_sig","gene_name","exon","rsid","dna_change","prot_change","transcript","zyg","inh","diseases_STR","cancer_groups"]]
    df_for_html.columns=["Classification","Gene","Exon/Intron","SNP ID","DNA change","Protein change","Transcript ID","Zygosity","Inheritance","Associated Disease(s)","Cancer Groups"]
    drug_response_df = df_for_html[df_for_html['Classification']=="drug response"]
    new_col = np.arange(1,(len(df_for_html) + 1))
    col_nbr_dict =  dict(zip(df_for_html["SNP ID"].tolist(),new_col))
    df_for_html.insert(0,"Number",new_col)
    drug_response_df.columns=["Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Response(s)","Cancer Groups"]
    df_for_html_2 = df_for_html[df_for_html['Classification']!="drug response"]
    df_for_html = df_for_html_2
    #print(df_for_html)
    #print(df_for_html["Classification"].tolist())
    #df_for_html['nbr'] = np.arange(len(df_for_html))
    #df_for_html = df_for_html[["nbr","clin_sig","gene_name","exon","rsid","dna_change","prot_change","transcript","zyg","inh","diseases_STR"]]
    #df_for_html.columns=["Number","Classification","Gene","Exon/\nIntron","SNP ID","DNA change","Protein\nchange","Transcript ID","Zygosity","Inheritance","Associated Disease(s)"]
    #print(list(ret_df.columns))
    df_for_freq = ret_df[["clin_sig","gene_name","transcript","dna_change","prot_change","chr_loc","freq","rsid"]]
    df_for_freq = df_for_freq[df_for_freq['clin_sig']!="drug response"]
    nbr_col_2 = [col_nbr_dict[rsid] for rsid in df_for_freq["rsid"].tolist()]
    df_for_freq = df_for_freq[["gene_name","transcript","dna_change","prot_change","chr_loc","freq","rsid"]]
    df_for_freq.columns=["Gene","Transcript","DNA change","Protein change","Genomic Location","Alternate Allele Fraction","dbSNP rsID"]
    #print(df_for_freq)
    #print(df_for_freq["rsid"].tolist())
    #print(nbr_col_2)
    list_of_str_to_remove = get_list_of_str_from_file(path_prefix + "dis_str_to_remove.txt")
    take_list_of_comments = "true"
    df_for_freq.insert(0,"Mutation\nNumber",nbr_col_2)
    comment_df = comment_df[comment_df['clin_sig']!="drug response"]
    comment_df = comment_df[comment_df['rsid'].notna()]
    nbr_col_3 = [col_nbr_dict[rsid] for rsid in comment_df["rsid"].tolist()]
    #print(comment_df["clin_sig"].tolist())
    comment_df.insert(0,"Number",nbr_col_3)
    if not(len(list(pubmed_df.index)) < 1):
        nbr_col_4 = [col_nbr_dict[rsid] for rsid in pubmed_df['rsid'].tolist()]
        pubmed_df.insert(0,"Number",nbr_col_4)
        pubmed_df = pubmed_df[["Number","name","PMID","title"]]
        pubmed_df.columns = ["Number","Name of Mutation","Pubmed ID","Title"]
    comments = comment_df.to_dict(orient='records')
    freqs = df_for_freq.to_dict(orient='records')
    pubmed = pubmed_df.to_dict(orient='records')
    #print(comments)
    #print(freqs)
    #print(pubmed)
    result_dict = df_for_html.to_dict(orient="records")
    dict_1 = {}
    ctr = 0
    dict_all = {}
    for gene_name_curr in list(df_for_html.index):
        curr_nbr_gene = str(df_for_html.loc[gene_name_curr,"Number"])
        gene_name_tmp = str(df_for_html.loc[gene_name_curr,"Gene"])
        gene_loc_tmp = str(df_for_html.loc[gene_name_curr,"DNA change"])
        gene_name = gene_name_tmp + " " + gene_loc_tmp
        new_dict = {}
        list_of_comments = []
        list_of_freqs = []
        list_of_pubmeds = []
        curr_row_as_dict = {}
        genes_found_by_cancer_type = {}
        for i in range(0,len(cancer_groups_for_filter)):
            genes_found_by_cancer_type[str(cancer_groups_for_filter[i])] = []
        if(cancer_groups_filter_active == "true"):
            curr_gene_cancer_groups = df_for_html.loc[gene_name_curr,"Cancer Groups"]
            for tmp in cancer_groups_for_filter:
                if(tmp in curr_gene_cancer_groups):
                    genes_found_by_cancer_type[tmp].append(gene_name)
            #if(found == "false"):
            #    continue
        if not(len(result_dict) < ctr):
            result_row = result_dict[ctr]
            ctr = ctr + 1
            diseases_for_curr_row = str(result_row["Associated Disease(s)"])
            disease_list_tmp = diseases_for_curr_row.split(";")
            s = []
            s_new = []
            for i in disease_list_tmp:
                #print(i)
                i_new = str("".join(i.split())).replace(" ","").replace("\t", "").lstrip()
                #print(i_new)
                #i_new = i_new.replace(" ","").replace("\t", "")
                if i_new.lower() not in s_new:
                    s.append(i)
                    s_new.append(i_new.lower())
                #print(s)
                #print(s_new)
            diseases_current = s
            disease_list_tmp = diseases_current
            #diseases_for_curr_row = str(result_row["Associated Disease(s)"])
            #disease_list_tmp = diseases_for_curr_row.split(";")
            result_row["Associated Disease(s)"] = disease_list_tmp
            if(gene_name in varsome_dict):
                result_row["Varsome"] = varsome_dict[gene_name]
            else:
                result_row["Varsome"] = {}
            curr_row_as_dict = result_row
        list_of_comments_dicts = []
        if(take_list_of_comments == "true"):
            comments_as_dict = []
        dict_of_existing_comments = {}
        for i in range(0,len(comments)):
            list_of_dicts = []
            if(take_list_of_comments == "true"):
                if not('Number' in comments[i]):
                    continue
                curr_nbr = str(comments[i]['Number'])
                if(curr_nbr == curr_nbr_gene):
                    del comments[i]["Number"]
                    del comments[i]["gene_name"]
                    del comments[i]["name"]
                    del comments[i]["rsid"]
                    print(comments[i])
                    #list_of_dicts = []
                    print(comments[i]['lab_name'])
                    print(len(comments[i]['lab_name']))
                    for j in range(0,len(comments[i]['lab_name'])):
                        print(j)
                        tmp_dict = {}
                        tmp_dict['clin_sig'] = comments[i]['clin_sig']
                        tmp_dict['lab_name'] = comments[i]['lab_name'][j]
                        tmp_dict['diseasename'] = comments[i]['diseasename'][j]
                        tmp_dict['comment'] = comments[i]['comment'][j]
                        tmp_dict['report_date'] = comments[i]['report_date'][j]
                        tmp_dict['labrank'] = comments[i]['labrank'][j]
                        dis_curr = str(comments[i]['diseasename'][j] + "")
                        dis_curr_new = str("".join(dis_curr.split())).replace(" ","").replace("\t", "").lstrip()
                        dis_curr_lower = dis_curr_new.lower()
                        if(dis_curr_lower in list_of_str_to_remove):
                            continue
                        if not(comments[i]['diseasename'][j] in dict_of_existing_comments):
                            dict_of_existing_comments[comments[i]['diseasename'][j]] = []
                        if not(comments[i]['lab_name'][j] in dict_of_existing_comments[comments[i]['diseasename'][j]]):
                            list_of_dicts.append(tmp_dict)
                        if not(comments[i]['diseasename'][j] in dict_of_existing_comments):
                            dict_of_existing_comments[comments[i]['diseasename'][j]] = [comments[i]['lab_name'][j]]
                        else:
                            dict_of_existing_comments[comments[i]['diseasename'][j]].append(comments[i]['lab_name'][j])
                        print(list_of_dicts)
                    #list_of_comments_dicts = list_of_dicts
                    print(list_of_comments_dicts)
                #list_of_comments_dicts = list_of_dicts
                list_of_comments_dicts = []
                list_of_diseasenames = [d['diseasename'] for d in list_of_dicts]
                list_of_diseasenames_unique = list(set(list_of_diseasenames))
                for diseasename in list_of_diseasenames_unique:
                    comments_with_this_disease = [d for d in list_of_dicts if (d['diseasename'] == diseasename)]
                    #sorted_list_curr = (sorted(comments_with_this_disease, key=lambda d: d['labrank']))
                    list_of_comments_dicts = list_of_comments_dicts + sorted(comments_with_this_disease, key=lambda d: d['labrank'])
                if(len(list_of_dicts) > 0):
                    comments_as_dict.append(list_of_comments_dicts)
                continue
            if not('Number' in comments[i]):
                continue
            curr_nbr = str(comments[i]['Number'])
            if(curr_nbr == curr_nbr_gene):
                del comments[i]["Number"]
                del comments[i]["gene_name"]
                del comments[i]["name"]
                del comments[i]["rsid"]
                tmp_dict = comments[i]
                list_of_comments.append(tmp_dict)
                comments_as_dict = tmp_dict
        for i in range(0,len(freqs)):
            if not('Mutation\nNumber' in freqs[i]):
                continue
            curr_nbr = str(freqs[i]['Mutation\nNumber'])
            if(curr_nbr == curr_nbr_gene):
                del freqs[i]["Mutation\nNumber"]
                tmp_dict = freqs[i]
                list_of_freqs.append(tmp_dict)
                frec_as_dict = tmp_dict
        for i in range(0,len(pubmed)):
            if not('Number' in pubmed[i]):
                continue
            curr_nbr = str(pubmed[i]['Number'])
            if(curr_nbr == curr_nbr_gene):
                del pubmed[i]["Number"]
                del pubmed[i]["Name of Mutation"]
                tmp_dict = pubmed[i]
                list_of_pubmeds.append(tmp_dict)
        new_dict['Gene Data'] = curr_row_as_dict
        new_dict['Comments'] = comments_as_dict
        new_dict['Frequencies'] = frec_as_dict
        new_dict['Publications'] = list_of_pubmeds
        dict_1[gene_name] = new_dict
        if(cancer_groups_filter_active == "true"):
            tmp = {}
            for group in cancer_groups_for_filter:
                if(gene_name in genes_found_by_cancer_type[group]):
                    tmp[gene_name] = new_dict
                    if not(group in dict_all):
                        dict_all[group] = [tmp]
                    else:
                        dict_all[group].append(tmp)
            dict_all_new = {}
            dict_all_new[customer_id] = dict_all
        dict_2 = {}
        dict_2[customer_id] = dict_1
    #results = df_for_html.to_dict(orient='records')
    if(cancer_groups_filter_active == "true"):
        results_json = json.dumps(dict_all_new)
    else:
        results_json = json.dumps(dict_2)
    results_json = results_json.replace("\n"," ").replace("NaN","\"NaN\"").replace("\"\"NaN\"\"","\"NaN\"")
    comments_json = json.dumps(comments)
    freqs_json = json.dumps(freqs)
    pubmed_json = json.dumps(pubmed)
    return(results_json,comments_json,freqs_json,pubmed_json)

def read_syndict(path):
    infile = open(path)
    htmlstr=infile.read()
    infile.close()
    lines= htmlstr.split("\n")
    syndict = {}
    for line in lines:
        lineSplit = line.split("\t")
        if(len(lineSplit) < 2):
            continue
        str1 = lineSplit[0]
        str2 = lineSplit[1]
        if not(str1 in syndict):
            syndict[str1] = []
        values = str2.split(",,")
        syndict[str1] = syndict[str1] + values
    return(syndict)


def read_cancer_groups(path):
    infile = open(path)
    htmlstr=infile.read()
    lines= htmlstr.split("\n")
    line_0 = lines[0]
    line0split = line_0.split("\t")
    cancer_groups = []
    genes_by_cancer_group = {}
    cancer_groups_by_gene = {}
    for i in range(1,len(line0split)):
        cancer_groups.append(line0split[i])
        genes_by_cancer_group[str(line0split[i])] = []
    ctr = 0
    for line in lines:
        if(ctr == 0):
            ctr = 1
            continue
        lineSplit = line.split("\t")
        gene_name = lineSplit[0]
        cancer_groups_by_gene[str(gene_name)] = []
        for i in range(1,len(lineSplit)):
            if(str(lineSplit[i]) == "TRUE"):
                genes_by_cancer_group[cancer_groups[i-1]].append(gene_name)
                cancer_groups_by_gene[gene_name].append(cancer_groups[i-1])
    #print(genes_by_cancer_group)
    #print(cancer_groups_by_gene)
    return([genes_by_cancer_group,cancer_groups_by_gene])

def read_varsome(rsid):
    #url = "https://api.varsome.com/lookup/rs2294081/1019?add-all-data=1&add-ACMG-annotation=1&annotation-mode=germline"
    #url = "https://api.varsome.com/lookup/rs10902121/1019?add-ACMG-annotation=1&add-all-data=1&annotation-mode=germline&format=api"
    url = "https://staging-api.varsome.com/lookup/chr17:7578443:A:T?add-AMP-annotation=1&format=api"
    # Making a get request
    user = "quirin@enigmagenomics.com" 
    passwd = "misterq.218" 
    # Make a request to the endpoint using the correct auth values 
    #auth_values = (user, passwd) 
    #response = requests.get(url, auth=auth_values) 
    #my_headers = {'Accept': 'application/json'}
    #auth99 = HTTPBasicAuth('apikey', '46tcqf#aw3hM0Qe!QFeZJeM&BZMpwkBl@oxbFfYC')
    #my_headers = {"Authorization" : '46tcqf#aw3hM0Qe!QFeZJeM&BZMpwkBl@oxbFfYC'}
    #response = requests.get(url, headers=my_headers_2)
    #response = requests.get(url,
    #                        auth = HTTPBasicAuth('quirin@enigmagenomics.com', 'misterq.218'), headers=my_headers)
    #response = urllib.request.urlopen(url)
    #cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/api-key -H 'content-type: application/x-www-form-urlencoded' -d apikey=b737cbc4-4006-44a1-b7c8-70e02102a7bd"
    cmd = "curl -X GET 'https://staging-api.varsome.com/lookup/" + rsid + "?add-AMP-annotation=1&format=api' -H  'Authorization: Token 46tcqf#aw3hM0Qe!QFeZJeM&BZMpwkBl@oxbFfYC' | iconv -f utf8 -t ascii//TRANSLIT//IGNORE"
    #cmd = "curl -X GET 'https://staging-api.varsome.com/lookup/rs28928907?add-AMP-annotation=1' -H  'Authorization: Token 46tcqf#aw3hM0Qe!QFeZJeM&BZMpwkBl@oxbFfYC'"
    stream = os.popen(cmd)
    webContent = stream.read()
    print(webContent)
    #return("none")
    #infile = open("varsome_test.txt")
    #webContent_tmp=infile.read()
    #infile.close()
    remove_if_term_not_found = "true"
    varsome_term_file = file_with_varsome_terms
    infile = open(varsome_term_file)
    varsome_str=infile.read()
    infile.close()
    list_of_varsome_terms = varsome_str.split("\n")
    #webContent = webContent_tmp
    #print(webContent)
    #gnomad_exomes = []
    #gnomad_genomes = []
    #gnomad_genomes_coverage = []
    ret = {}
    for curr_element in list_of_varsome_terms:
        if(len(curr_element) < 2):
            if curr_element in list_of_varsome_terms:
                list_of_varsome_terms.remove(curr_element)
    for term_curr in list_of_varsome_terms:
        ret[term_curr] = []
    #ret["gnomad_exomes"] = []
    #ret["gnomad_genomes"] = []
    #ret["gnomad_genomes_coverage"] = []
    if not("<b>Vary:</b> <span class=\"lit\">Accept</span>" in webContent):
        return(ret)
    tmp_1 = webContent.split("<b>Vary:</b> <span class=\"lit\">Accept</span>")[1]
    tmp_2 = tmp_1.split("</span>")[1]
    #print(tmp_2)
    tmp_3 = tmp_2.split("</pre>")[0]
    tmp4 = tmp_3.strip("[").rstrip("]")
    #json_all = json.loads(tmp_3)[0]
    #print(tmp4)
    json_all_tmp = json.loads(tmp_3)
    #print("json ")
    #print(json_all_tmp)
    json_all = []
    if not(isinstance(json_all_tmp, dict)):
        if(isinstance(json_all_tmp, list)):
            json_all = json_all_tmp[0]
            #print("json_all is list")
    else:
        json_all = json_all_tmp
        #print("json_all is dict")
    #print("keys:")
    #for key in json_all_tmp:
    #    print(key)
    #        #print(json_all[key]["gnomad_exomes"])
    #    #print(json_all[key])
    for term_curr in list_of_varsome_terms:
        if(term_curr in json_all):
            curr_entry_for_this_varsome_term = json_all[term_curr]
            ret[term_curr] = curr_entry_for_this_varsome_term
        elif(remove_if_term_not_found == "true"):
            ret.pop(term_curr, None)
    #if("gnomad_exomes" in json_all):
    #    gnomad_exomes = json_all["gnomad_exomes"]
    #if("gnomad_genomes" in json_all):
    #    gnomad_genomes = json_all["gnomad_genomes"]
    #if("gnomad_genomes_coverage" in json_all):
    #    gnomad_genomes_coverage = json_all["gnomad_genomes_coverage"]
    #ret["gnomad_exomes"] = gnomad_exomes
    #ret["gnomad_genomes"] = gnomad_genomes
    #ret["gnomad_genomes_coverage"] = gnomad_genomes_coverage
    #print(json_all)
    return(ret)


def get_list_of_str_from_file(path):
    infile = open(path)
    htmlstr=infile.read()
    lines= htmlstr.split("\n")
    list_to_remove = []
    ctr = 0
    list_to_remove_ctr_list = []
    for line in lines:
        line_new = str("".join(line.split())).replace(" ","").replace("\t", "").lstrip()
        line_lower = line_new.lower()
        list_to_remove.append(line_lower)
        list_to_remove.append(line)
    return(list_to_remove)

def remove_str_from_list(path,list_of_dis):
    infile = open(path)
    htmlstr=infile.read()
    lines= htmlstr.split("\n")
    list_to_remove = []
    ctr = 0
    list_to_remove_ctr_list = []
    for line in lines:
        line_new = str("".join(line.split())).replace(" ","").replace("\t", "").lstrip()
        line_lower = line_new.lower()
        list_to_remove.append(line_lower)
    #list_of_dis_new = [str("".join(dis.split())).replace(" ","").replace("\t", "").lstrip() for dis in list_of_dis]
    #list_of_dis_lower = [dis.lower() for dis in list_of_dis_new]
    #dis_new = [dis for dis in list_of_dis_lower if not(dis in list_to_remove)]
    dis_for_ret = []
    for i in range(0,len(list_of_dis)):
        dis_curr = list_of_dis[i]
        dis_curr_old = str(dis_curr + "")
        dis_curr_new = str("".join(dis_curr.split())).replace(" ","").replace("\t", "").lstrip()
        dis_curr_lower = dis_curr_new.lower()
        if not(dis_curr_lower in list_to_remove):
            dis_for_ret.append(dis_curr_old)
    return(dis_for_ret)

def remove_str_from_dis_and_sig(path,dis_and_sig):
    infile = open(path)
    htmlstr=infile.read()
    lines= htmlstr.split("\n")
    list_to_remove = []
    for line in lines:
        line_new = str("".join(line.split())).replace(" ","").replace("\t", "").lstrip()
        line_lower = line_new.lower()
        list_to_remove.append(line_lower)
    dis_and_sig_as_list = list(dis_and_sig)
    ret = []
    for i in dis_and_sig_as_list:
        #print(str(i))
        dis_old = str(i + "")
        dis_curr = str("".join(str(i).split(";")[0].split())).replace(" ","").replace("\t", "").lstrip()
        dis_curr_lower = dis_curr.lower()
        if not(dis_curr_lower in list_to_remove):
            ret.append(dis_old)
    return(ret)
            
            
def remove_synonyms(list_of_diseases):
    syn = [get_synonyms(i,"TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas") for i in list_of_diseases]
    syn_exists = [i for i in syn if len(i) > 0]
    if(len(syn_exists) < 1):
        return(list_of_diseases)
    for i in range(0,len(list_of_diseases)):
        for j in range(0,len(list_of_diseases)):
            if(i == j or (list_of_diseases[i] == "") or (list_of_diseases[j] == "")):
                continue
            else:
                if (list_of_diseases[j].lower() in syn[i]):
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    list_of_diseases[j] = ""
                elif(list_of_diseases[i].lower() in syn[j]):
                    #print("removed" + str(list_of_diseases[i])+ "(dupl. of" + str(list_of_diseases[j]) + ")")
                    list_of_diseases[i] = ""
                elif(list_of_diseases[i].lower() == list_of_diseases[j].lower()):
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    list_of_diseases[j] = ""
                else:
                    continue
    list_of_diseases_2 = list_of_diseases.copy()
    list_of_diseases = [i for i in list_of_diseases_2 if not(i == "")]
    ret = list(dict.fromkeys(list_of_diseases))
    #for i in range(0,len(list_of_diseases)):
    #    if(list_of_diseases[i] == ""):
    #        list_of_diseases_2.pop(i)
    return(ret)

def remove_synonyms_2(dis_and_sig):
    syn = [get_synonyms((i.split(";")[0]),"TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas") for i in dis_and_sig]
    syn_exists = [i for i in syn if len(i) > 0]
    if(len(syn_exists) < 1):
        return(dis_and_sig)
    dis_and_sig_list = list(dis_and_sig)
    list_of_diseases = [i.split(";")[0] for i in dis_and_sig]
    list_of_significances = [i.split(";")[1] for i in dis_and_sig]
    if(len(dis_and_sig_list) < 1):
        return(dict(zip([],[])))
    #print(list(dis_and_sig.values()))
    for i in range(0,len(list_of_diseases)):
        for j in range(0,len(list_of_diseases)):
            #print(list(dis_and_sig.values())[i])
            if(pd.isnull(list(dis_and_sig.values())[i]) or pd.isnull(list(dis_and_sig.values())[j])):
                continue
            if(i == j or (list_of_diseases[i] == "") or (list_of_diseases[j] == "")):
                continue
            else:
                if (list_of_diseases[j].lower() in syn[i] or list_of_diseases[i].lower() in syn[j] or list_of_diseases[i].lower() == list_of_diseases[j].lower()):
                    if("N/A" in list(dis_and_sig.values())[i]):
                        #print("removed element :" + str(list_of_diseases[i]))
                        list_of_diseases[i] = ""
                    elif("N/A" in list(dis_and_sig.values())[j]):
                        #print("removed element :" + str(list_of_diseases[j]))
                        list_of_diseases[j] = ""
                    elif("pathogenic" in list_of_significances[i] and not("pathogenic" in list_of_significances[j])):
                        #print("removed element :" + str(list_of_diseases[j]))
                        list_of_diseases[j] = ""
                    elif("pathogenic" in list_of_significances[j] and not("pathogenic" in list_of_significances[i])):
                        #print("removed element :" + str(list_of_diseases[i]))
                        list_of_diseases[i] = ""
                    
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    #list_of_diseases[j] = ""
                else:
                    continue
    list_of_diseases_2 = list_of_diseases.copy()
    ret_arr_keys = [dis_and_sig_list[i] for i in range(0,len(dis_and_sig_list)) if not(list_of_diseases_2[i] == "")]
    ret_arr_values = [dis_and_sig[dis_and_sig_list[i]] for i in range(0,len(dis_and_sig_list)) if not(list_of_diseases_2[i] == "")]
    ret = dict(zip(ret_arr_keys,ret_arr_values))
    #list_of_diseases = [i for i in list_of_diseases_2 if not(i == "")]
    #ret = list(dict.fromkeys(list_of_diseases))
    #for i in range(0,len(list_of_diseases)):
    #    if(list_of_diseases[i] == ""):
    #        list_of_diseases_2.pop(i)
    return(ret)

def remove_synonyms_2_new(dis_and_sig,syndict):
    #syn = [get_synonyms((i.split(";")[0]),"TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas") for i in dis_and_sig]
    #syn_exists = [i for i in syn if len(i) > 0]
    #if(len(syn_exists) < 1):
    #    return(dis_and_sig)
    syn = []
    list_of_diseases_new = [i.split(";")[0] for i in dis_and_sig]
    for i in range(0, len(list_of_diseases_new)):
        curr_list = list_of_diseases_new[i]
        [syn_curr,syndict_2] = get_synonyms_new(curr_list,"TGT-1324012-CDXW7ErncSwjrLVeYAKog5b2oMJUZtkPV3DXMYNer9l1kNs5sQ-cas",syndict)
        syndict = syndict_2.copy()
        syn.append(syn_curr)
        #print(curr_list)
        #print("syn ")
        #print(syn_curr)
    dis_and_sig_list = list(dis_and_sig)
    list_of_diseases = [i.split(";")[0] for i in dis_and_sig]
    list_of_significances = [i.split(";")[1] for i in dis_and_sig]
    if(len(dis_and_sig_list) < 1):
        return(dict(zip([],[])))
    #print(list(dis_and_sig.values()))
    for i in range(0,len(list_of_diseases)):
        for j in range(0,len(list_of_diseases)):
            #print(list(dis_and_sig.values())[i])
            if(pd.isnull(list(dis_and_sig.values())[i]) or pd.isnull(list(dis_and_sig.values())[j])):
                continue
            if(i == j or (list_of_diseases[i] == "") or (list_of_diseases[j] == "")):
                continue
            else:
                if (list_of_diseases[j].lower() in syn[i] or list_of_diseases[i].lower() in syn[j] or list_of_diseases[i].lower() == list_of_diseases[j].lower()):
                    if("N/A" in list(dis_and_sig.values())[i]):
                        #print("removed element :" + str(list_of_diseases[i]))
                        list_of_diseases[i] = ""
                    elif("N/A" in list(dis_and_sig.values())[j]):
                        #print("removed element :" + str(list_of_diseases[j]))
                        list_of_diseases[j] = ""
                    elif("pathogenic" in list_of_significances[i] and not("pathogenic" in list_of_significances[j])):
                        #print("removed element :" + str(list_of_diseases[j]))
                        list_of_diseases[j] = ""
                    elif("pathogenic" in list_of_significances[j] and not("pathogenic" in list_of_significances[i])):
                        #print("removed element :" + str(list_of_diseases[i]))
                        list_of_diseases[i] = ""
                    
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    #list_of_diseases[j] = ""
                else:
                    continue
    list_of_diseases_2 = list_of_diseases.copy()
    ret_arr_keys = [dis_and_sig_list[i] for i in range(0,len(dis_and_sig_list)) if not(list_of_diseases_2[i] == "")]
    ret_arr_values = [dis_and_sig[dis_and_sig_list[i]] for i in range(0,len(dis_and_sig_list)) if not(list_of_diseases_2[i] == "")]
    ret = dict(zip(ret_arr_keys,ret_arr_values))
    #list_of_diseases = [i for i in list_of_diseases_2 if not(i == "")]
    #ret = list(dict.fromkeys(list_of_diseases))
    #for i in range(0,len(list_of_diseases)):
    #    if(list_of_diseases[i] == ""):
    #        list_of_diseases_2.pop(i)
    return(ret,syndict_2)
    
def remove_synonyms_3(dis_and_sig):
    syn = [get_synonyms((i.split(";")[0]),"TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas") for i in dis_and_sig]
    syn_exists = [i for i in syn if len(i) > 0]
    if(len(syn_exists) < 1):
        return(dis_and_sig)
    dis_and_sig_list = list(dis_and_sig)
    list_of_diseases = [i.split(";")[0] for i in dis_and_sig]
    list_of_significances = [i.split(";")[1] for i in dis_and_sig]
    if(len(dis_and_sig_list) < 1):
        return(dict(zip([],[])))
    #print(list(dis_and_sig.values()))
    for i in range(0,len(list_of_diseases)):
        for j in range(0,len(list_of_diseases)):
            #print(list(dis_and_sig.values())[i])
            if(pd.isnull(list(dis_and_sig.values())[i]) or pd.isnull(list(dis_and_sig.values())[j])):
                continue
            if(i == j or (list_of_diseases[i] == "") or (list_of_diseases[j] == "")):
                continue
            else:
                if ((list_of_diseases[j].lower() in syn[i] or list_of_diseases[i].lower() in syn[j] or list_of_diseases[i].lower() == list_of_diseases[j].lower()) and (list_of_significances[i] == list_of_significances[j])):
                    if("N/A" in list(dis_and_sig.values())[i]):
                        #print("removed element :" + str(list_of_diseases[i]))
                        list_of_diseases[i] = ""
                    elif("N/A" in list(dis_and_sig.values())[j]):
                        #print("removed element :" + str(list_of_diseases[j]))
                        list_of_diseases[j] = ""
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    #list_of_diseases[j] = ""
                else:
                    continue
    list_of_diseases_2 = list_of_diseases.copy()
    ret_arr_keys = [dis_and_sig_list[i] for i in range(0,len(dis_and_sig_list)) if not(list_of_diseases_2[i] == "")]
    ret_arr_values = [dis_and_sig[dis_and_sig_list[i]] for i in range(0,len(dis_and_sig_list)) if not(list_of_diseases_2[i] == "")]
    ret = dict(zip(ret_arr_keys,ret_arr_values))
    #list_of_diseases = [i for i in list_of_diseases_2 if not(i == "")]
    #ret = list(dict.fromkeys(list_of_diseases))
    #for i in range(0,len(list_of_diseases)):
    #    if(list_of_diseases[i] == ""):
    #        list_of_diseases_2.pop(i)
    return(ret)
    
        

def remove_synonyms_new(list_of_diseases,syndict):
    #[syn,syndicts_2] = [get_synonyms_new(i,"TGT-1324012-CDXW7ErncSwjrLVeYAKog5b2oMJUZtkPV3DXMYNer9l1kNs5sQ-cas",syndict) for i in list_of_diseases]
    syn = []
    for i in range(0, len(list_of_diseases)):
        curr_list = list_of_diseases[i]
        [syn_curr,syndict_2] = get_synonyms_new(curr_list,"TGT-1324012-CDXW7ErncSwjrLVeYAKog5b2oMJUZtkPV3DXMYNer9l1kNs5sQ-cas",syndict)
        syndict = syndict_2.copy()
        syn.append(syn_curr)
        #print(curr_list)
        #print("syn ")
        #print(syn_curr)
    #syndict_2 = syndicts_2[len(syndicts_2)-1]
    #syn = [get_synonyms_new(i,"TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas",syndict) for i in list_of_diseases]
    for i in range(0,len(list_of_diseases)):
        for j in range(0,len(list_of_diseases)):
            if(i == j or (list_of_diseases[i] == "") or (list_of_diseases[j] == "")):
                continue
            else:
                if (list_of_diseases[j].lower() in syn[i]):
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    list_of_diseases[j] = ""
                elif(list_of_diseases[i].lower() in syn[j]):
                    #print("removed" + str(list_of_diseases[i])+ "(dupl. of" + str(list_of_diseases[j]) + ")")
                    list_of_diseases[i] = ""
                elif(list_of_diseases[i].lower() == list_of_diseases[j].lower()):
                    #print("removed" + str(list_of_diseases[j]) + "(dupl. of" + str(list_of_diseases[i]) + ")")
                    list_of_diseases[j] = ""
                else:
                    continue
    list_of_diseases_2 = list_of_diseases.copy()
    list_of_diseases = [i for i in list_of_diseases_2 if not(i == "")]
    ret_2 = list(dict.fromkeys(list_of_diseases))
    #for i in range(0,len(list_of_diseases)):
    #    if(list_of_diseases[i] == ""):
    #        list_of_diseases_2.pop(i)
    return(ret_2,syndict_2)


def read_cancer_groups_OLD(path):
    infile = open(path)
    htmlstr=infile.read()
    lines= htmlstr.split("\n")
    line_0 = lines[0]
    line0split = line_0.split("\t")
    cancer_groups = []
    genes_by_cancer_group = {}
    cancer_groups_by_gene = {}
    for i in range(1,len(line0split)):
        cancer_groups.append(line0split[i])
        genes_by_cancer_group[str(line0split[i])] = []
    ctr = 0
    for line in lines:
        if(ctr == 0):
            ctr = 1
            continue
        lineSplit = line.split("\t")
        gene_name = lineSplit[0]
        cancer_groups_by_gene[str(gene_name)] = []
        for i in range(1,len(lineSplit)):
            if(str(lineSplit[i]) == "TRUE"):
                genes_by_cancer_group[cancer_groups[i-1]].append(gene_name)
                cancer_groups_by_gene[gene_name].append(cancer_groups[i-1])
    #print(genes_by_cancer_group)
    #print(cancer_groups_by_gene)


def calculate_polygenic(prsice_path,disease,gwas_res_path,plink_path,bg_plink_dataset):
    return()
def extract_vep_data_OLD(dataframe,disease_filter,start,end):
    df = dataframe
    #print(df)
    pubmed_list = pd.DataFrame(columns=['name','PMID','title'])
    [gen_dis_dict,dis_gen_dict] = read_gene_associations(assoc_table_path)
    [gene_dis_dict_2,dis_gene_dict_2] = read_gene_associations_2(gene_disease_groups_path)
    publ_ctr = 0
    hgnc_dict = {}
    inh_dict = {}
    gene_info_dict = {}
    snp_dict = {}
    ret_df = pd.DataFrame(columns=['rsid','clin_sig','consequence','variant_class','dna_change','prot_change','exon','transcript','chr_loc','hgvs','gene_name','hgnc','freq','max_af','max_pop','vaf','gq','ad1','ad2','ad','zyg','inh','clin_sig_list_STR','clin_sig_list','clin_sig_ct','diseases','diseases_STR','disease_groups','comments','comments_STR','dis_and_sig'])
    syndict_all = read_syndict("/scripts/syndict_temp.txt")
    #for i in range(0,len(df.index)):
    #for i in range(1400,1600):
    #for i in range(2000,3000):
    #print(dis_gene_dict_2[disease])
    #print(gene_dis_dict_2["EPHA2"])
    #print([str(i) for i in gene_dis_dict_2])
    #print([i for i in gene_dis_dict_2 if disease_filter in gene_dis_dict_2[i]])
    #if(use_all_lines == "true"):
    #    start_line_def = 0
    #    stop_line_def = 19000
    start_line_def = int(start)
    stop_line_def = int(end)
    for i in range(start_line_def,min(len(df.index),stop_line_def)):
        rsid_tmp = df.loc[i,"Existing_variation"]
        temp1 = rsid_tmp.split(",")
        rsids = [i for i in temp1 if("rs" in i)]
        rsids_unique = [i for i in temp1 if(("rs" in i) and not(i in df.index))]
        if(any(rsid_curr in rsids for rsid_curr in list(df.index.values))):
            if(len(rsids_unique) > 0):
                rsid = rsids_unique[0]
            else:
                continue
        if(len(rsids) > 0):
            rsid = rsids[0]
            variant_class = df.loc[i,"VARIANT_CLASS"]
        else:
            continue
        gene_name = df.loc[i,"SYMBOL"]
        if ((not(disease_filter in dis_gene_dict_2)) and (not(disease_filter=="all"))):
            #print("Disease group does not exist! Returning empty dataframes.")
            return([pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()])
        if not(disease_filter=="all"):
            if not(gene_name in gene_dis_dict_2):
                #print(str(gene_name) + "not in gene_dis_dict")
                continue
            #else:
            #    print(gene_dis_dict_2[gene_name])
            #print(str(disease))
            elif(not(str(disease_filter) in ",".join(gene_dis_dict_2[gene_name]))):
                #print(str(disease_filter))
                #print(gene_dis_dict_2[gene_name])
                continue
        #print(i)
        hgnc_id = str(df.loc[i,"HGNC_ID"])
        if not(rsid in hgnc_dict):
            hgnc_dict[rsid] = hgnc_id
        prot_change = df.loc[i,"HGVSp"]
        if(prot_change == "-" or not("ENSP" in prot_change) or (len(prot_change.split(":")) < 2)):
            continue
        prot_change = prot_change.split(":")[1]
        prot_change = prot_change.replace("%3D","=")
        consequence = df.loc[i,"Consequence"]
        consequence = consequence.split(",")[0]
        #exon = str(df.loc[i,"EXON"])
        #if(exon == "-"):
        #    exon = str(df.loc[i,"INTRON"])
        #exon = str(exon.split("/")[0])
        exon = "1"
        chr_loc = df.loc[i,"Location"]
        if not("Chr" in chr_loc or "chr" in chr_loc):
            chr_loc = "Chr"+str(chr_loc)
        af = df.loc[i,"gnomAD_AF"]
        max_af = df.loc[i,"MAX_AF"]
        if(str(max_af) == "-"):
            max_af = "N/A"
        elif(not(max_af.isdigit()) and str(max_af).replace('.','',1).isdigit()):
            max_af = float(max_af)
            max_af = str("{:.2%}".format(max_af))
        elif(str(max_af) == "0" or str(max_af) == "1"):
            max_af = float(max_af)
            max_af = str("{:.2%}".format(max_af))
        else:
            max_af = "N/A"
        
        max_pop = df.loc[i,"MAX_AF_POPS"]
        max_pop = max_pop.replace("gnomAD_","")
        if not(max_pop == "-"):
            max_pop_arr = max_pop.split(",")
            max_pop_arr_2 = list(dict.fromkeys(max_pop_arr))
            max_pop_arr_3 = []
            for pop in max_pop_arr_2:
                if(pop in pop_dict):
                    max_pop_arr_3.append(pop_dict[pop])
                else:
                    max_pop_arr_3.append(pop)
            max_pop = ",".join(max_pop_arr_3)
        else:
            max_pop = "N/A"
        vaf = df.loc[i,"Sample_1.VAF"]
        if not(str(vaf) == "NA"):
            if(str(vaf).replace('.','',1).isdigit()):
                vaf = float(vaf)
                vaf = str("{:.2%}".format(vaf))
        else:
            vaf = "N/A"
        gq = df.loc[i,"Sample_1.GQ"]
        if(str(gq).isdigit()):
            gq = float(gq)
        ad = df.loc[i,"Sample_1.AD"]
        #print(str(ad))
        #print(str(ad).split(","))
        ad1 = ""
        ad2 = ""
        if(str(ad).split(",")[0].isdigit() and (len(str(ad).split(",")) > 1)):
            ad = df.loc[i,"Sample_1.AD"]
            ad1 = float(str(ad).split(",")[0])
            ad2 = float(str(ad).split(",")[1])
            ad = int(round(ad1 + ad2))
        else:
            ad = "NA"
            #ad = str(ad1) + "/" + str(ad2)
        zyg_tmp = df.loc[i,"Sample_1.GT"]
        #print(zyg_tmp)
        if(zyg_tmp.split("/")[0] == "." or (len(zyg_tmp.split("/")) < 1)):
            zyg = "---"
            continue
        elif(str(zyg_tmp.split("/")[0]) == str(zyg_tmp.split("/")[1])):
            zyg = "homozygous"
        else:
            zyg = "heterozygous"
        #print(zyg)
        if(str(af) == "-"):
            af = 999.0
        else:
            af = float(str(af))
            af = str("{:.2%}".format(af))
        #if(hgnc_id in inh_dict):
        #    inh = inh_dict[hgnc_id]
        #else:
        #    inh = get_inheritance_mode_2(hgnc_id)
        #    inh_dict[hgnc_id] = inh
        if(gene_name in inh_dict):
            inh = inh_dict[gene_name]
        else:
            inh = get_inheritance_mode_3(gene_name,"/scripts/tableExport.csv","","false")
            inh_dict[gene_name] = inh
        #inh = get_inheritance_mode(gene_name)
        #print(inh)
        #print(prot_change)
        clin_sig_conf_tmp = df.loc[i,"ClinVar_updated_2019Dec_CLNSIGCONF"]
        clin_sig_list_STR = ""
        clin_sig_list = []
        clin_sig_ct = []
        if not(str(clin_sig_conf_tmp) == "-"):
            sig_list = clin_sig_conf_tmp.split("%3B")
            if(len(sig_list) > 1):
                for n in sig_list:
                    clin_sig_list_STR = clin_sig_list_STR + ";" + str(n.split("(")[0])
                    clin_sig_list.append(str(n.split("(")[0]))
                    clin_sig_ct.append(str(n.split("(")[1].split(")")[0]))
            else:
                clin_sig_list_STR = cln_sig_conf_tmp.split("(")[0]
        clin_sig_2 = ""
        if(len(clin_sig_list) > 1):
            index_main_sig = clin_sig_ct.index(max(clin_sig_ct))
            if(index_main_sig < len(clin_sig_list)):
                clin_sig_2 = clin_sig_list[index_main_sig]
        clin_sig = df.loc[i,"CLIN_SIG"]
        clin_sig.replace("_"," ")
        has_clinvar = "false"
        if("," in clin_sig):
            if(len(str(df.loc[i,"ClinVar_updated_2019Dec_CLNSIG"])) > 2):
                clin_sig = df.loc[i,"ClinVar_updated_2019Dec_CLNSIG"]
                clin_sig = clin_sig.replace("_"," ")
        if(clin_sig == "-"):
            polyphen = df.loc[i,"PolyPhen"]
            if(")" in polyphen):
                #print("polyphen")
                #print(polyphen)
                polyphen_sig = polyphen.split("(")[0]
                polyphen_val = float(polyphen.split("(")[1].split(")")[0])
                clin_sig = polyphen_sig.replace("probably_damaging","likely pathogenic")
                clin_sig = clin_sig.replace("_"," ")
                #print(str(polyphen_val))
                #print(clin_sig)
                if not (clin_sig == "benign" or clin_sig == "likely pathogenic" or (polyphen_val > polyphen_threshold)):
                    continue
            else:
                continue
        else:
            has_clinvar = "true"
        if(clin_sig_2 == ""):
            clin_sig_2 = clin_sig
        gene_ensembl_id = df.loc[i,"Gene"]
        if("ENSG" in gene_ensembl_id and not( gene_name in gene_info_dict)):
            [gene_entrez_id,summary] = get_gene_info(gene_ensembl_id)
            gene_info_dict[gene_name] = summary
        disease = df.loc[i,"ClinVar_updated_2019Dec_CLNDN"]
        #if(disease == "-"):
        #    continue
        if not(disease == "-"):
            disease = disease.replace("_"," ")
            diseases_curr_rs = disease.split("|")
            #print(i)
            #print(diseases_curr_rs)
            diseases_curr_rs_tmp = [j for j in diseases_curr_rs if not(("not specified" in j) or ("not provided" in j) or ("not_specified" in j) or ("not_provided" in j))]
            if(len(diseases_curr_rs_tmp) > 1):
                #diseases_curr_rs = remove_synonyms(diseases_curr_rs_tmp)
                [diseases_curr_rs,syndict_2] = remove_synonyms_new(diseases_curr_rs_tmp,syndict_all)
                syndict = syndict_2.copy()
                syndict_all.update(syndict_2)
            else:
                diseases_curr_rs = diseases_curr_rs_tmp
        else:
            diseases_curr_rs = []
        #print(diseases_curr_rs)
        #for j in range(0,arrlen):
        #    if(("not specified" in diseases_curr_rs[j]) or ("not provided" in diseases_curr_rs[j]) or ("not_specified" in diseases_curr_rs[j]) or ("not_provided" in diseases_curr_rs[j])):
        #        diseases_curr_rs.pop(j)
        #        j = j - 1
        #        arrlen = arrlen - 1
        #print(disease)
        if(has_clinvar == "true"):
            hgvs_list = get_hgvs_list(rsid)
        else:
            hgvs_list = []
        gene_name_dict[rsid] = gene_name
        prot_change_dict[rsid] = prot_change
        snp_data_tmp = {}
        hgvsc = df.loc[i,"HGVSc"]
        if(":" in hgvsc):
            transcript = hgvsc.split(":")[0]
            dna_change = hgvsc.split(":")[1]
        else:
            transcript = "-"
            dna_change = "-"
        snp_data_tmp = {}
        if (rsid in list(df.index.values)):
            continue
        if not(gene_name_dict[rsid] in snp_dict):
            snp_dict[gene_name_dict[rsid]] = []
        snp_dict[gene_name_dict[rsid]].append(rsid)
        snp_data_tmp["clin_sig"] = clin_sig
        snp_data_tmp["dna_change"] = dna_change
        snp_data_tmp["prot_change"] = prot_change
        ret_df.loc[rsid,"clin_sig"] = clin_sig
        ret_df.loc[rsid,"dna_change"] = dna_change
        ret_df.loc[rsid,"prot_change"] = prot_change
        ret_df.loc[rsid,"freq"] = af
        ret_df.loc[rsid,"max_af"] = max_af
        ret_df.loc[rsid,"max_pop"] = max_pop
        ret_df.loc[rsid,"vaf"] = vaf
        #print(gq)
        #print(rsid)
        #print(inh)
        ret_df.loc[rsid,"gq"] = gq
        ret_df.loc[rsid,"ad1"] = ad1
        ret_df.loc[rsid,"ad2"] = ad2
        ret_df.loc[rsid,"ad"] = ad
        ret_df.loc[rsid,"hgnc_id"] = hgnc_id
        ret_df.loc[rsid,"inh"] = inh
        ret_df.loc[rsid,"exon"] = exon
        ret_df.loc[rsid,"chr_loc"] = chr_loc
        ret_df.loc[rsid,"consequence"] = consequence
        ret_df.loc[rsid,"transcript"] = transcript
        ret_df.loc[rsid,"zyg"] = zyg
        ret_df.loc[rsid,"clin_sig_list_STR"] = clin_sig_list_STR
        ret_df.loc[rsid,"clin_sig_list"] = clin_sig_list
        ret_df.loc[rsid,"clin_sig_ct"] = clin_sig_ct
        ret_df.loc[rsid,"rsid"] = rsid
        ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
        [hgvs_curr,link] = get_correct_hgvs(hgvs_list,gene_name,prot_change)
        #ret = extract(link)
        if(len(hgvs_list) < 1 or (hgvs_curr == "")):
            if("benign" in clin_sig and remove_benign == "true" and not("drug response" in clin_sig)):
                continue
            hgvs_curr = "-"
            #print("SNP ID:" + str(rsid) + "; Gene Name: " + gene_name + "\n")
            #print("Associated disease(s):\n")
            #print(disease)
            #if(rsid in dict_by_snpid):
            #    for dis in diseases_curr_rs:
            #        dict_by_snpid[rsid][disease] = "N/A"
            #else:
            #    for dis in diseases_curr_rs:
            #        dict_by_snpid[rsid] = {}
            #        dict_by_snpid[rsid][dis] = "N/A"
            #for dis in diseases_curr_rs:
            #    if(dis in ret_all):
            #        ret_all[dis][rsid]="N/A"
            #    else:
            #        ret_all[dis] = {}
            #        ret_all[dis][rsid]="N/A"
            if(rsid in dict_by_snpid):
                for dis in diseases_curr_rs:
                    dict_by_snpid[rsid][str(dis + ";" + clin_sig)] = "N/A"
            else:
                for dis in diseases_curr_rs:
                    dict_by_snpid[rsid] = {}
                    dict_by_snpid[rsid][str(dis + ";" + clin_sig)] = "N/A"
            for dis in diseases_curr_rs:
                if(str(dis + ";" + clin_sig) in ret_all):
                    ret_all[str(dis + ";" + clin_sig)][rsid]="N/A"
                else:
                    ret_all[str(dis + ";" + clin_sig)] = {}
                    ret_all[str(dis + ";" + clin_sig)][rsid]="N/A"
            if(rsid not in dict_by_snpid):
                continue
            ret_df.loc[rsid,"dis_and_sig"] = [k for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"diseases"] = [k.split(";")[0] for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"diseases_STR"] = ";".join([k.split(";")[0] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments"] = [dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]]
            if not(gene_name in gene_dis_dict_2):
                ret_df.loc[rsid,"disease_groups"] = "N/A"
            else:
                ret_df.loc[rsid,"disease_groups"] = ";".join(gene_dis_dict_2[gene_name])
            #ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
            #ret_df.loc[rsid,"chr_loc"] = chr_loc
            #ret_df.loc[rsid,"rsid"] = rsid
            #ret_df.loc[rsid,"exon"] = exon
            #ret_df.loc[rsid,"inh"] = inh
            #ret_df.loc[rsid,"freq"] = af
            #print([dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments_STR"] = ";".join([str(dict_by_snpid[rsid][k]) for k in dict_by_snpid[rsid]])
            if(len(diseases_curr_rs) > 1):
                snp_data_tmp["hgvs"] = "-"
                ret_df.loc[rsid,"hgvs"] = "-"
                snp_data[rsid]=snp_data_tmp
            else:
                continue
        else:
            #print("SNP ID:" + str(rsid) + "; Gene Name: " + gene_name + "\n")
            #print("Associated disease(s):\n")
            #[hgvs_curr,link] = get_correct_hgvs(hgvs_list,gene_name,prot_change)
            ret = extract(link)
            dis_and_sig_from_vep_keys = [(str(dis) + ";" + str(clin_sig)) for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            dis_and_sig_from_vep_values = ["N/A" for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            if("benign" in clin_sig and remove_benign=="true" and not("drug response" in clin_sig)):
                dis_and_sig_from_vep_keys = []
                dis_and_sig_from_vep_values = []
            dis_and_sig_curr = dict(zip(dis_and_sig_from_vep_keys,dis_and_sig_from_vep_values))
            #dis_and_sig_curr_rs = [(str(dis) + ";" + str(clin_sig)) for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            #ret2 = remove_synonyms_2[ret]
            #print("lists: ")
            #print(dis_and_sig_from_vep_keys)
            #print(list(ret))
            for key in list(ret):
                if(("not specified" in key) or ("not provided" in key) or ("not_specified" in key) or ("not_provided" in key)):
                    del ret[key]
                elif("benign" in key and remove_benign=="true" and not("drug response" in clin_sig)):
                    del ret[key]
                else:
                    dis_and_sig_curr[key] = ret[key]
                    #dis_and_sig_curr_rs.append(key)
            #print(i)
            [dis_and_sig_curr_rs,syndict_2] = remove_synonyms_2_new(dis_and_sig_curr,syndict_all)
            syndict = syndict_2.copy()
            syndict_all.update(syndict_2)
            #print(dis_and_sig_curr_rs)
            #diseases_curr_rs_tmp = [j for j in diseases_curr_rs if not(("not specified" in j) or ("not provided" in j) or ("not_specified" in j) or ("not_provided" in j))]
            #print("d.split(;)[0]")
            #print([d.split(";")[0] for d in dis_and_sig_curr_rs])
            #dis_no_syn = remove_synonyms([d.split(";")[0] for d in dis_and_sig_curr_rs])
            #print(dis_no_syn)
            #dis_and_sig_curr_rs = [dis for dis in dis_and_sig_curr_rs if (dis.split(";")[0] in dis_no_syn)]
            if(len(ret) < 1):
                continue
            for dis in dis_and_sig_curr_rs:
                if(rsid in dict_by_snpid):
                    dict_by_snpid[rsid][dis] = dis_and_sig_curr_rs[dis]
                    #dict_by_snpid[rsid][dis] = ret[dis]
                else:
                    dict_by_snpid[rsid] = {}
                    dict_by_snpid[rsid][dis] = dis_and_sig_curr_rs[dis]
                    #dict_by_snpid[rsid][dis] = ret[dis]
                if(dis in ret_all):
                    ret_all[dis][rsid]= dis_and_sig_curr_rs[dis]
                    #ret_all[dis][rsid]=ret[dis]
                else:
                    ret_all[dis] = {}
                    ret_all[dis][rsid]= dis_and_sig_curr_rs[dis]
                    #ret_all[dis][rsid]=ret[dis]
                if(dis not in dis_and_sig_from_vep_keys):
                    if(rsid not in dict_by_snpid):
                        continue
            #print(dict_by_snpid[rsid])
            #snp_data_tmp[rsid]=dict_by_snpid[rsid]
            snp_data_tmp["hgvs"] = hgvs_curr
            ret_df.loc[rsid,"hgvs"] = hgvs_curr
            snp_data_tmp["gene_name"] = gene_name
            ret_df.loc[rsid,"gene_name"] = gene_name
            #snp_data_tmp["diseases"] = ";".join([key for key in ret])
            snp_data[rsid]=snp_data_tmp
            #print(ret_df)
            ret_df.loc[rsid,"diseases"] = [k.split(";")[0] for k in dict_by_snpid[rsid]]
            #disease_with_inh = ""
            #for disease_name_curr in dict_by_snpid[rsid]:
            #    inh = get_inheritance_mode_3(gene_name,"VEP_test/tableExport.csv",disease_name_curr,"false")
            ret_df.loc[rsid,"diseases_STR"] = ";".join([k.split(";")[0] for k in dict_by_snpid[rsid]])
            if not(gene_name in gene_dis_dict_2):
                ret_df.loc[rsid,"disease_groups"] = "N/A"
            else:
                ret_df.loc[rsid,"disease_groups"] = ";".join(gene_dis_dict_2[gene_name])
            ret_df.loc[rsid,"comments"] = [dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]]
            #print([dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments_STR"] = ";".join([str(dict_by_snpid[rsid][k]) for k in dict_by_snpid[rsid]])
            #ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
            ret_df.loc[rsid,"dis_and_sig"] = [k for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"consequence"] = consequence
            #ret_df.loc[rsid,"transcript"] = transcript
            #ret_df.loc[rsid,"inh"] = inh
            #ret_df.loc[rsid,"freq"] = af
            #ret_df.loc[rsid,"exon"] = exon
            #ret_df.loc[rsid,"chr_loc"] = chr_loc
            #ret_df.loc[rsid,"zyg"] = zyg
            ret_df.loc[rsid,"clin_sig_list_STR"] = clin_sig_list_STR
            ret_df.loc[rsid,"clin_sig_list"] = clin_sig_list
            ret_df.loc[rsid,"clin_sig_ct"] = clin_sig_ct
            [pubmed_dict,pmid_STR_tmp,pmid_STR] = get_publications(rsid)
            if not(pmid_STR == "-9999"):
                curr_pubmed_ctr = 0
                for key in pubmed_dict:
                    if not(curr_pubmed_ctr > pubmed_df_max_nr):
                        name = str(ret_df.loc[rsid,"gene_name"]) + " " + str(ret_df.loc[rsid,"dna_change"]) + " (" + str(ret_df.loc[rsid,"prot_change"]) + ")"
                        pubmed_list.loc[publ_ctr,"name"] = name
                        pubmed_list.loc[publ_ctr,"rsid"] = rsid
                        pubmed_list.loc[publ_ctr,"PMID"] = key
                        pubmed_list.loc[publ_ctr,"title"] = pubmed_dict[key]
                        publ_ctr = publ_ctr + 1
                        curr_pubmed_ctr = curr_pubmed_ctr + 1
            if(len(ret_df.loc[rsid,"dis_and_sig"]) < 2):
                if(pd.isnull(ret_df.loc[rsid,"dis_and_sig"])):
                    ret_df.loc[rsid,"dis_and_sig"] = []
            ret_df.loc[rsid,"rsid"] = rsid
    #for key in ret_all:
    #    print("Disease Name: " + str(key))
    #    for key2 in ret_all[key]:
    #        print("SNP ID:" + key2 + "; Gene Name: " + gene_name_dict[key2])
    #        print(ret_all[key][key2])
    #print(dict_by_snpid)
    #for snp_curr in snp_data:
        #print(snp_data[snp_curr])
        #print("SNP ID: " + str(snp_curr))
        #print(snp_data[snp_curr])
        #print(gene_name_dict[snp_curr])
        #print(dict_by_snpid[snp_curr])
        #print("Classification:\t" + snp_data[snp_curr]["clin_sig"] + "\tGene:\t" + gene_name_dict[snp_curr] + "\tDNA change:\t" + snp_data[snp_curr]["dna_change"] + "\tProtein change:\t" + snp_data[snp_curr]["prot_change"] + "\tAssociated disease(s):\t" + ";".join([i.split(";")[0] for i in dict_by_snpid[snp_curr]]))
        #for dis_and_sig in dict_by_snpid[snp_curr]:
        #    if(len(str(dis_and_sig).split(";")) > 1):
        #           print("RSID:\t" + snp_curr + "\t - \t" + str(dis_and_sig).split(";")[0] + ";\t" + str(dis_and_sig).split(";")[1] +":")
        #           if not(str(dict_by_snpid[snp_curr][dis_and_sig]) == "nan"):
        #               print("Description:\t" + str(dict_by_snpid[snp_curr][dis_and_sig]))
    comment_df = pd.DataFrame(columns=['clin_sig','name','comment','rsid'])
    ctr = 0
    ctr2 = 0
    comment_htmls = {}
    comment_and_disease_dict = {}
    list_of_described_genes = []
    for i in list(ret_df.index):
        #print(str(ret_df.loc[i,"rsid"]))
        #print(str(ret_df.loc[i,"dis_and_sig"]))
        if(isinstance(ret_df.loc[i,"dis_and_sig"], list)):
            if(len(ret_df.loc[i,"dis_and_sig"]) < 2):
                if(str(ret_df.loc[i,"rsid"]) == "nan" or pd.isnull(ret_df.loc[i,"rsid"]) or str(ret_df.loc[i,"rsid"]) == "NaN" or (str(ret_df.loc[i,"diseases"]) == "nan" or pd.isnull(ret_df.loc[i,"diseases"]) or str(ret_df.loc[i,"diseases"]) == "NaN")):
                    #print("dropped" + str(ret_df.loc[i,"rsid"]))
                    ret_df = ret_df.drop([i])
                    continue
        else:
            if(str(ret_df.loc[i,"rsid"]) == "nan" or pd.isnull(ret_df.loc[i,"rsid"]) or str(ret_df.loc[i,"rsid"]) == "NaN" or (str(ret_df.loc[i,"dis_and_sig"]) == "nan" or pd.isnull(ret_df.loc[i,"dis_and_sig"]) or str(ret_df.loc[i,"dis_and_sig"]) == "NaN") or (str(ret_df.loc[i,"diseases"]) == "nan" or pd.isnull(ret_df.loc[i,"diseases"]) or str(ret_df.loc[i,"diseases"]) == "NaN")):
                ret_df = ret_df.drop([i])
                continue
        if not(str(ret_df.loc[i,"comments_STR"]) == "nan"):
            if((len(ret_df.loc[i,"comments_STR"]) > 2) and not (str(ret_df.loc[i,"comments_STR"]) == "nan") and not (str(ret_df.loc[i,"comments_STR"])=="N/A")):
                    for k in range(0, len(list(ret_df.loc[i,"dis_and_sig"]))):
                        if not (str(list(ret_df.loc[i,"comments"])[k]) == "nan" or str(list(ret_df.loc[i,"comments"])[k]) == "N/A"):
                            name = str(ret_df.loc[i,"gene_name"]) + " " + str(ret_df.loc[i,"dna_change"]) + " (" + str(ret_df.loc[i,"prot_change"]) + ") - "
                            name = name + str(list(ret_df.loc[i,"diseases"])[k])
                            comment_df.loc[name,"clin_sig"] = str(ret_df.loc[i,"clin_sig"])
                            comment_df.loc[name,"rsid"] = str(ret_df.loc[i,"rsid"])
                            comment_df.loc[name,"comment"] = str(list(ret_df.loc[i,"comments"])[k])
                            comment_df.loc[name,"gene_name"] = str(ret_df.loc[i,"gene_name"])
        #if not(str(ret_df.loc[i,"gene_name"]) in gen_dis_dict):
        #    gen_dis_dict[str(ret_df.loc[i,"gene_name"])] = list(ret_df.loc[i,"diseases"])
        #else:
        #    gen_dis_dict[str(ret_df.loc[i,"gene_name"])].append(list(ret_df.loc[i,"diseases"]))
    clin_sig_all_values = ret_df["clin_sig"].tolist()
    clin_sig_all_values_STR = "".join(clin_sig_all_values)
    is_pathogenic = "false"
    if("pathogenic" in clin_sig_all_values_STR):
        is_pathogenic = "true"
    disease_list = []
    disease_dict_for_ret = {}
    for i in list(ret_df.index):
        dis_tmp_1 = ret_df.loc[i,"diseases"]
        disease_list = dis_tmp_1
        if(ret_df.loc[i,"gene_name"] in gen_dis_dict):
            dis_tmp_2 = gen_dis_dict[ret_df.loc[i,"gene_name"]]
            disease_list = dis_tmp_1 + dis_tmp_2
        [new_disease_list,syndict_2] = remove_synonyms_new(disease_list,syndict_all)
        disease_dict_for_ret[ret_df.loc[i,"gene_name"]] = new_disease_list
        syndict = syndict_2.copy()
        syndict_all.update(syndict_2)
    gene_info_ret = {}
    for gene in gene_info_dict:
        if(gene in disease_dict_for_ret):
            #gene_info_ret[gene] = [gene_info_dict[gene],disease_dict_for_ret[gene]]
            gene_info_ret[gene] = [gene_info_dict[gene],disease_dict_for_ret[gene],snp_dict[gene]]
    #for i in list(comment_df.index):
    #    gene = comment_df.loc[i,"gene_name"]
    #    if (not(gene in list_of_described_genes) and (comment_df.loc[i,"gene_name"] in gene_info_dict)):
    #        comment_htmls[ctr2] = gene_info_dict[comment_df.loc[i,"gene_name"]]
    #        if not(i in comment_and_disease_dict):
    #            comment_and_disease_dict[i] = ctr2
    #        else:
    #            comment_and_disease_dict[i].append(ctr2)
    #        ctr2 = ctr2 + 1
    #    if not(i in comment_and_disease_dict):
    #        comment_and_disease_dict[i] = ctr2
    #    else:
    #        comment_and_disease_dict[i].append(ctr2)
    #    comment_htmls[ctr2] = comment_df.loc[i,"comment"]
    return([ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic])
    #return([ret_df,comment_df,pubmed_list])
    #return[ret_all,dict_by_snpid,snp_data]
    #print(df.head())



def extract_vep_data(dataframe,disease_filter,start,end,sample_id):
    #df = pd.read_csv(path,sep='\t')
    df = dataframe
    #print(df)
    pubmed_list = pd.DataFrame(columns=['name','PMID','title'])
    [gen_dis_dict,dis_gen_dict] = read_gene_associations(assoc_table_path)
    [gene_dis_dict_2,dis_gene_dict_2] = read_gene_associations_2(gene_disease_groups_path)
    publ_ctr = 0
    hgnc_dict = {}
    inh_dict = {}
    gene_info_dict = {}
    snp_dict = {}
    ret_df = pd.DataFrame(columns=['rsid','clin_sig','consequence','variant_class','dna_change','prot_change','exon','transcript','chr_loc','hgvs','gene_name','hgnc','freq','max_af','max_pop','vaf','gq','ad1','ad2','ad','zyg','inh','clin_sig_list_STR','clin_sig_list','clin_sig_ct','diseases','diseases_STR','disease_groups','comments','comments_STR','dis_and_sig','cancer_groups'])
    samplename = sample_id
    syndict_all = read_syndict(str(path_prefix + "syndict_temp.txt"))
    cancer_groups_tmp = read_cancer_groups(str(path_prefix +"genes_and_cancer_groups.txt"))
    [genes_by_cancer_group,cancer_groups_by_gene] = cancer_groups_tmp
    #for i in range(0,len(df.index)):
    #for i in range(1400,1600):
    #for i in range(2000,3000):
    #print(dis_gene_dict_2[disease])
    #print(gene_dis_dict_2["EPHA2"])
    #print([str(i) for i in gene_dis_dict_2])
    #print([i for i in gene_dis_dict_2 if disease_filter in gene_dis_dict_2[i]])
    if(use_all_lines == "true"):
        start_line_def = 0
        stop_line_def = 20000
    if(automatic_stop_line_def == "true"):
        if(len(df.index) > 6700):
            start_line_def = 6600
            stop_line_def = 6700
        elif(len(df.index) > 100):
            start_line_def = 0
            stop_line_def = 100
        else:
            start_line_def = start
            stop_line_def = end
    else:
        start_line_def = start
        stop_line_def = end
    varsome_dict = {}
    clinvar_version = "ClinVar_updated_2021Aug"
    if("ClinVar_updated_2019Dec_CLNDN" in df):
        clinvar_version = "ClinVar_updated_2019Dec"
    for i in range(start_line_def,min(len(df.index),stop_line_def)):
        rsid_tmp = df.loc[i,"Existing_variation"]
        temp1 = rsid_tmp.split(",")
        rsids = [i for i in temp1 if("rs" in i)]
        rsids_unique = [i for i in temp1 if(("rs" in i) and not(i in df.index))]
        if(any(rsid_curr in rsids for rsid_curr in list(df.index.values))):
            if(len(rsids_unique) > 0):
                rsid = rsids_unique[0]
            else:
                continue
        if(len(rsids) > 0):
            rsid = rsids[0]
            variant_class = df.loc[i,"VARIANT_CLASS"]
        else:
            continue
        gene_name = df.loc[i,"SYMBOL"]
        if ((not(disease_filter in dis_gene_dict_2)) and (not(disease_filter=="all"))):
            #print("Disease group does not exist! Returning empty dataframes.")
            return([pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()])
        if not(disease_filter=="all"):
            if not(gene_name in gene_dis_dict_2):
                #print(str(gene_name) + "not in gene_dis_dict")
                continue
            #else:
            #    print(gene_dis_dict_2[gene_name])
            #print(str(disease))
            elif(not(str(disease_filter) in ",".join(gene_dis_dict_2[gene_name]))):
                #print(str(disease_filter))
                #print(gene_dis_dict_2[gene_name])
                continue
        #print(i)
        hgnc_id = str(df.loc[i,"HGNC_ID"])
        if not(rsid in hgnc_dict):
            hgnc_dict[rsid] = hgnc_id
        prot_change = df.loc[i,"HGVSp"]
        if(prot_change == "-" or not("ENSP" in prot_change) or (len(prot_change.split(":")) < 2)):
            continue
        prot_change = prot_change.split(":")[1]
        prot_change = prot_change.replace("%3D","=")
        consequence = df.loc[i,"Consequence"]
        consequence = consequence.split(",")[0]
        #exon = str(df.loc[i,"EXON"])
        #if(exon == "-"):
        #    exon = str(df.loc[i,"INTRON"])
        #exon = str(exon.split("/")[0])
        exon = "1"
        chr_loc = df.loc[i,"Location"]
        if not("Chr" in chr_loc or "chr" in chr_loc):
            chr_loc = "Chr"+str(chr_loc)
        af = df.loc[i,"gnomAD_AF"]
        max_af = df.loc[i,"MAX_AF"]
        if(str(max_af) == "-"):
            max_af = "N/A"
        elif(not(str(max_af).isdigit()) and str(max_af).replace('.','',1).isdigit()):
            max_af = float(str(max_af))
            max_af = str("{:.2%}".format(max_af))
        elif(str(max_af) == "0" or str(max_af) == "1"):
            max_af = float(str(max_af))
            max_af = str("{:.2%}".format(max_af))
        else:
            max_af = "N/A"
        
        max_pop = df.loc[i,"MAX_AF_POPS"]
        max_pop = max_pop.replace("gnomAD_","")
        if not(max_pop == "-"):
            max_pop_arr = max_pop.split(",")
            max_pop_arr_2 = list(dict.fromkeys(max_pop_arr))
            max_pop_arr_3 = []
            for pop in max_pop_arr_2:
                if(pop in pop_dict):
                    max_pop_arr_3.append(pop_dict[pop])
                else:
                    max_pop_arr_3.append(pop)
            max_pop = ",".join(max_pop_arr_3)
        else:
            max_pop = "N/A"
        vaf = df.loc[i, samplename +".VAF"]
        if not(str(vaf) == "NA"):
            if(str(max_af).replace('.','',1).isdigit()):
                vaf = float(str(vaf))
                vaf = str("{:.2%}".format(vaf))
        else:
            vaf = "N/A"
        #if("Sample_1.GQ" in df):
        gq = df.loc[i,samplename +".GQ"]
        #elif("GQ" in df):
        #    vaf = df.loc[i,"GQ"]
        #else:
        #    gq = 0.0
        #gq = df.loc[i,"Sample_1.GQ"]
        if(str(gq).isdigit()):
            gq = float(str(gq))
        ad = df.loc[i,samplename +".AD"]
        #print(str(ad))
        #print(str(ad).split(","))
        ad1 = ""
        ad2 = ""
        if(str(ad).split(",")[0].isdigit() and (len(str(ad).split(",")) > 1)):
            ad = df.loc[i,samplename +".AD"]
            ad1 = float(str(ad).split(",")[0])
            ad2 = float(str(ad).split(",")[1])
            ad = int(round(ad1 + ad2))
        else:
            ad = "NA"
            #ad = str(ad1) + "/" + str(ad2)
        zyg_tmp = df.loc[i,samplename +".GT"]
        #print(zyg_tmp)
        if(zyg_tmp.split("/")[0] == "." or (len(zyg_tmp.split("/")) < 1)):
            zyg = "---"
            continue
        elif(str(zyg_tmp.split("/")[0]) == str(zyg_tmp.split("/")[1])):
            zyg = "homozygous"
        else:
            zyg = "heterozygous"
        #print(zyg)
        if(str(af) == "-"):
            af = 999.0
        else:
            af = float(str(af))
            af = str("{:.2%}".format(af))
        #if(hgnc_id in inh_dict):
        #    inh = inh_dict[hgnc_id]
        #else:
        #    inh = get_inheritance_mode_2(hgnc_id)
        #    inh_dict[hgnc_id] = inh
        #if(hgnc_id in inh_dict):
        #    inh = inh_dict[hgnc_id]
        #else:
        #    inh = get_inheritance_mode_2(hgnc_id)
        #    inh_dict[hgnc_id] = inh
        if(gene_name in inh_dict):
            inh = inh_dict[gene_name]
        else:
            inh = get_inheritance_mode_3(gene_name,str(path_prefix + "tableExport.csv"),"","false")
            inh_dict[gene_name] = inh
        #inh = get_inheritance_mode(gene_name)
        #print(inh)
        #print(prot_change)
        clin_sig_conf_tmp = df.loc[i,clinvar_version + "_CLNSIGCONF"]
        clin_sig_list_STR = ""
        clin_sig_list = []
        clin_sig_ct = []
        if not(str(clin_sig_conf_tmp) == "-"):
            sig_list = clin_sig_conf_tmp.split("%3B")
            if(len(sig_list) > 1):
                for n in sig_list:
                    clin_sig_list_STR = clin_sig_list_STR + ";" + str(n.split("(")[0])
                    clin_sig_list.append(str(n.split("(")[0]))
                    clin_sig_ct.append(str(n.split("(")[1].split(")")[0]))
            else:
                clin_sig_list_STR = clin_sig_conf_tmp.split("(")[0]
        clin_sig_2 = ""
        if(len(clin_sig_list) > 1):
            index_main_sig = clin_sig_ct.index(max(clin_sig_ct))
            if(index_main_sig < len(clin_sig_list)):
                clin_sig_2 = clin_sig_list[index_main_sig]
        clin_sig = df.loc[i,"CLIN_SIG"]
        clin_sig.replace("_"," ")
        has_clinvar = "false"
        if("," in clin_sig):
            if(len(str(df.loc[i,clinvar_version + "_CLNSIG"])) > 2):
                clin_sig = df.loc[i,clinvar_version + "_CLNSIG"]
                clin_sig = clin_sig.replace("_"," ")
        if(clin_sig == "-"):
            polyphen = df.loc[i,"PolyPhen"]
            if(")" in polyphen):
                #print("polyphen")
                #print(polyphen)
                polyphen_sig = polyphen.split("(")[0]
                polyphen_val = float(str(polyphen.split("(")[1].split(")")[0]))
                clin_sig = polyphen_sig.replace("probably_damaging","likely pathogenic")
                clin_sig = clin_sig.replace("_"," ")
                #print(str(polyphen_val))
                #print(clin_sig)
                if not (clin_sig == "benign" or clin_sig == "likely pathogenic" or (polyphen_val > polyphen_threshold)):
                    #print("No clinical significance, polyphen value given")
                    continue
            else:
                #print("No clinical significance")
                continue
        else:
            has_clinvar = "true"
        if(clin_sig_2 == ""):
            clin_sig_2 = clin_sig
        gene_ensembl_id = df.loc[i,"Gene"]
        if("ENSG" in gene_ensembl_id and not( gene_name in gene_info_dict)):
            [gene_entrez_id,summary] = get_gene_info(gene_ensembl_id)
            gene_info_dict[gene_name] = summary
        disease = df.loc[i,clinvar_version + "_CLNDN"]
        #if(disease == "-"):
        #    continue
        if not(disease == "-"):
            disease = disease.replace("_"," ")
            diseases_curr_rs = disease.split("|")
            #print(diseases_curr_rs)
            #print(i)
            #print(diseases_curr_rs)
            diseases_curr_rs_tmp = [j for j in diseases_curr_rs if not(("not specified" in j) or ("not provided" in j) or ("not_specified" in j) or ("not_provided" in j))]
            if(len(diseases_curr_rs_tmp) > 1):
                diseases_curr_rs = remove_synonyms(diseases_curr_rs_tmp)
            else:
                diseases_curr_rs = diseases_curr_rs_tmp
        else:
            diseases_curr_rs = []
        #print(diseases_curr_rs)
        #for j in range(0,arrlen):
        #    if(("not specified" in diseases_curr_rs[j]) or ("not provided" in diseases_curr_rs[j]) or ("not_specified" in diseases_curr_rs[j]) or ("not_provided" in diseases_curr_rs[j])):
        #        diseases_curr_rs.pop(j)
        #        j = j - 1
        #        arrlen = arrlen - 1
        #print(disease)
        if(has_clinvar == "true"):
            hgvs_list = get_hgvs_list(rsid)
        else:
            hgvs_list = []
        gene_name_dict[rsid] = gene_name
        prot_change_dict[rsid] = prot_change
        snp_data_tmp = {}
        hgvsc = df.loc[i,"HGVSc"]
        if(":" in hgvsc):
            transcript = hgvsc.split(":")[0]
            dna_change = hgvsc.split(":")[1]
        else:
            transcript = "-"
            dna_change = "-"
        snp_data_tmp = {}
        gene_name_copy = "" + str(gene_name)
        gene_name_uppercase = gene_name_copy.upper()
        cancer_groups = []
        if(gene_name_uppercase in cancer_groups_by_gene):
            cancer_groups = cancer_groups_by_gene[gene_name_uppercase]
        if (rsid in list(df.index.values)):
            continue
        if not(gene_name_dict[rsid] in snp_dict):
            snp_dict[gene_name_dict[rsid]] = []
        snp_dict[gene_name_dict[rsid]].append(rsid)
        snp_data_tmp["clin_sig"] = clin_sig
        snp_data_tmp["dna_change"] = dna_change
        snp_data_tmp["prot_change"] = prot_change
        ret_df.loc[rsid,"clin_sig"] = clin_sig
        ret_df.loc[rsid,"dna_change"] = dna_change
        ret_df.loc[rsid,"prot_change"] = prot_change
        ret_df.loc[rsid,"freq"] = af
        ret_df.loc[rsid,"max_af"] = max_af
        ret_df.loc[rsid,"max_pop"] = max_pop
        ret_df.loc[rsid,"vaf"] = vaf
        ret_df.loc[rsid,"cancer_groups"] = cancer_groups
        #print(gq)
        #print(rsid)
        #print(inh)
        ret_df.loc[rsid,"gq"] = gq
        ret_df.loc[rsid,"ad1"] = ad1
        ret_df.loc[rsid,"ad2"] = ad2
        ret_df.loc[rsid,"ad"] = ad
        ret_df.loc[rsid,"hgnc_id"] = hgnc_id
        ret_df.loc[rsid,"inh"] = inh
        ret_df.loc[rsid,"exon"] = exon
        ret_df.loc[rsid,"chr_loc"] = chr_loc
        ret_df.loc[rsid,"consequence"] = consequence
        ret_df.loc[rsid,"transcript"] = transcript
        ret_df.loc[rsid,"zyg"] = zyg
        ret_df.loc[rsid,"clin_sig_list_STR"] = clin_sig_list_STR
        ret_df.loc[rsid,"clin_sig_list"] = clin_sig_list
        ret_df.loc[rsid,"clin_sig_ct"] = clin_sig_ct
        ret_df.loc[rsid,"rsid"] = rsid
        ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
        [hgvs_curr,link] = get_correct_hgvs(hgvs_list,gene_name,prot_change)
        #ret = extract(link)
        if(len(hgvs_list) < 1 or (hgvs_curr == "")):
            if("benign" in clin_sig and remove_benign == "true" and not("drug response" in clin_sig)):
                continue
            hgvs_curr = "-"
            #print("SNP ID:" + str(rsid) + "; Gene Name: " + gene_name + "\n")
            #print("Associated disease(s):\n")
            #print(disease)
            #if(rsid in dict_by_snpid):
            #    for dis in diseases_curr_rs:
            #        dict_by_snpid[rsid][disease] = "N/A"
            #else:
            #    for dis in diseases_curr_rs:
            #        dict_by_snpid[rsid] = {}
            #        dict_by_snpid[rsid][dis] = "N/A"
            #for dis in diseases_curr_rs:
            #    if(dis in ret_all):
            #        ret_all[dis][rsid]="N/A"
            #    else:
            #        ret_all[dis] = {}
            #        ret_all[dis][rsid]="N/A"
            if(rsid in dict_by_snpid):
                for dis in diseases_curr_rs:
                    dict_by_snpid[rsid][str(dis + ";" + clin_sig)] = "N/A"
            else:
                for dis in diseases_curr_rs:
                    dict_by_snpid[rsid] = {}
                    dict_by_snpid[rsid][str(dis + ";" + clin_sig)] = "N/A"
            for dis in diseases_curr_rs:
                if(str(dis + ";" + clin_sig) in ret_all):
                    ret_all[str(dis + ";" + clin_sig)][rsid]="N/A"
                else:
                    ret_all[str(dis + ";" + clin_sig)] = {}
                    ret_all[str(dis + ";" + clin_sig)][rsid]="N/A"
            if(rsid not in dict_by_snpid):
                continue
            ret_df.loc[rsid,"dis_and_sig"] = [k for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"diseases"] = [k.split(";")[0] for k in dict_by_snpid[rsid]]
            curr_dna_change = ret_df.loc[rsid,"dna_change"]
            gene_name_tmp = ret_df.loc[rsid,"gene_name"]
            gene_name_and_mutation = gene_name_tmp + " " + curr_dna_change
            varsome_dict[gene_name_and_mutation] = read_varsome(rsid)
            ret_df.loc[rsid,"diseases_STR"] = ";".join([k.split(";")[0] for k in dict_by_snpid[rsid]])
            curr_dis = list(ret_df.loc[rsid,"diseases"])
            curr_dis_new = remove_str_from_list((path_prefix + "dis_str_to_remove.txt"),curr_dis)
            curr_dis_tmp_2 = str(ret_df.loc[rsid,"diseases_STR"]).split(";")
            curr_dis_tmp_2_new = remove_str_from_list((path_prefix + "dis_str_to_remove.txt"),curr_dis_tmp_2)
            dis_and_sig_curr_temp = ret_df.loc[rsid,"dis_and_sig"]
            dis_and_sig_new = remove_str_from_dis_and_sig((path_prefix + "dis_str_to_remove.txt"),dis_and_sig_curr_temp)
            ret_df.loc[rsid,"dis_and_sig"] = dis_and_sig_new
            if(len(curr_dis_new) == 0):
                ret_df = ret_df.drop([rsid])
                continue
            else:
                ret_df.loc[rsid,"diseases"] = curr_dis_new
                ret_df.loc[rsid,"diseases_STR"] = ";".join(curr_dis_tmp_2_new)
            ret_df.loc[rsid,"comments"] = [dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]]
            if not(gene_name in gene_dis_dict_2):
                ret_df.loc[rsid,"disease_groups"] = "N/A"
            else:
                ret_df.loc[rsid,"disease_groups"] = ";".join(gene_dis_dict_2[gene_name])
            #ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
            #ret_df.loc[rsid,"chr_loc"] = chr_loc
            #ret_df.loc[rsid,"rsid"] = rsid
            #ret_df.loc[rsid,"exon"] = exon
            #ret_df.loc[rsid,"inh"] = inh
            #ret_df.loc[rsid,"freq"] = af
            #print([dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments_STR"] = ";".join([str(dict_by_snpid[rsid][k]) for k in dict_by_snpid[rsid]])
            if(len(diseases_curr_rs) > 1):
                snp_data_tmp["hgvs"] = "-"
                ret_df.loc[rsid,"hgvs"] = "-"
                snp_data[rsid]=snp_data_tmp
            else:
                continue
        else:
            #print("SNP ID:" + str(rsid) + "; Gene Name: " + gene_name + "\n")
            #print("Associated disease(s):\n")
            #[hgvs_curr,link] = get_correct_hgvs(hgvs_list,gene_name,prot_change)
            ret = extract(link)
            dis_and_sig_from_vep_keys = [(str(dis) + ";" + str(clin_sig)) for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            dis_and_sig_from_vep_values = ["N/A" for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            if("benign" in clin_sig and remove_benign=="true" and not("drug response" in clin_sig)):
                dis_and_sig_from_vep_keys = []
                dis_and_sig_from_vep_values = []
            dis_and_sig_curr = dict(zip(dis_and_sig_from_vep_keys,dis_and_sig_from_vep_values))
            #dis_and_sig_curr_rs = [(str(dis) + ";" + str(clin_sig)) for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            #ret2 = remove_synonyms_2[ret]
            #print("lists: ")
            #print(dis_and_sig_from_vep_keys)
            #print(list(ret))
            for key in list(ret):
                if(("not specified" in key) or ("not provided" in key) or ("not_specified" in key) or ("not_provided" in key)):
                    del ret[key]
                elif("benign" in key and remove_benign=="true" and not("drug response" in clin_sig)):
                    del ret[key]
                else:
                    dis_and_sig_curr[key] = ret[key]
                    #dis_and_sig_curr_rs.append(key)
            #print(i)
            [dis_and_sig_curr_rs,syndict_2] = remove_synonyms_2_new(dis_and_sig_curr,syndict_all)
            syndict = syndict_2.copy()
            syndict_all.update(syndict_2)
            #print(dis_and_sig_curr_rs)
            #diseases_curr_rs_tmp = [j for j in diseases_curr_rs if not(("not specified" in j) or ("not provided" in j) or ("not_specified" in j) or ("not_provided" in j))]
            #print("d.split(;)[0]")
            #print([d.split(";")[0] for d in dis_and_sig_curr_rs])
            #dis_no_syn = remove_synonyms([d.split(";")[0] for d in dis_and_sig_curr_rs])
            #print(dis_no_syn)
            #dis_and_sig_curr_rs = [dis for dis in dis_and_sig_curr_rs if (dis.split(";")[0] in dis_no_syn)]
            if(len(ret) < 1):
                continue
            for dis in dis_and_sig_curr_rs:
                if(rsid in dict_by_snpid):
                    dict_by_snpid[rsid][dis] = dis_and_sig_curr_rs[dis]
                    #dict_by_snpid[rsid][dis] = ret[dis]
                else:
                    dict_by_snpid[rsid] = {}
                    dict_by_snpid[rsid][dis] = dis_and_sig_curr_rs[dis]
                    #dict_by_snpid[rsid][dis] = ret[dis]
                if(dis in ret_all):
                    ret_all[dis][rsid]= dis_and_sig_curr_rs[dis]
                    #ret_all[dis][rsid]=ret[dis]
                else:
                    ret_all[dis] = {}
                    ret_all[dis][rsid]= dis_and_sig_curr_rs[dis]
                    #ret_all[dis][rsid]=ret[dis]
                if(dis not in dis_and_sig_from_vep_keys):
                    if(rsid not in dict_by_snpid):
                        continue
            #print(dict_by_snpid[rsid])
            #snp_data_tmp[rsid]=dict_by_snpid[rsid]
            snp_data_tmp["hgvs"] = hgvs_curr
            ret_df.loc[rsid,"hgvs"] = hgvs_curr
            snp_data_tmp["gene_name"] = gene_name
            ret_df.loc[rsid,"gene_name"] = gene_name
            #snp_data_tmp["diseases"] = ";".join([key for key in ret])
            snp_data[rsid]=snp_data_tmp
            #print(ret_df)
            ret_df.loc[rsid,"diseases"] = [k.split(";")[0] for k in dict_by_snpid[rsid]]
            #disease_with_inh = ""
            #for disease_name_curr in dict_by_snpid[rsid]:
            #    inh = get_inheritance_mode_3(gene_name,"VEP_test/tableExport.csv",disease_name_curr,"false")
            ret_df.loc[rsid,"diseases_STR"] = ";".join([k.split(";")[0] for k in dict_by_snpid[rsid]])
            if not(gene_name in gene_dis_dict_2):
                ret_df.loc[rsid,"disease_groups"] = "N/A"
            else:
                ret_df.loc[rsid,"disease_groups"] = ";".join(gene_dis_dict_2[gene_name])
            ret_df.loc[rsid,"comments"] = [dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]]
            #print([dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments_STR"] = ";".join([str(dict_by_snpid[rsid][k]) for k in dict_by_snpid[rsid]])
            #ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
            curr_dna_change = ret_df.loc[rsid,"dna_change"]
            gene_name_tmp = ret_df.loc[rsid,"gene_name"]
            gene_name_and_mutation = gene_name_tmp + " " + curr_dna_change
            varsome_dict[gene_name_and_mutation] = read_varsome(rsid)
            ret_df.loc[rsid,"dis_and_sig"] = [k for k in dict_by_snpid[rsid]]
            #print(ret_df.loc[rsid,"dis_and_sig"])
            curr_dis = list(ret_df.loc[rsid,"diseases"])
            curr_dis_new = remove_str_from_list((path_prefix + "dis_str_to_remove.txt"),curr_dis)
            curr_dis_tmp_2 = str(ret_df.loc[rsid,"diseases_STR"]).split(";")
            curr_dis_tmp_2_new = remove_str_from_list((path_prefix + "dis_str_to_remove.txt"),curr_dis_tmp_2)     
            dis_and_sig_curr_temp = ret_df.loc[rsid,"dis_and_sig"]
            dis_and_sig_new = remove_str_from_dis_and_sig((path_prefix + "dis_str_to_remove.txt"),dis_and_sig_curr_temp)
            ret_df.loc[rsid,"dis_and_sig"] = dis_and_sig_new
            if(len(curr_dis_new) == 0):
                ret_df = ret_df.drop([rsid])
                continue
            else:
                ret_df.loc[rsid,"diseases"] = curr_dis_new
                ret_df.loc[rsid,"diseases_STR"] = ";".join(curr_dis_tmp_2_new)
            ret_df.loc[rsid,"consequence"] = consequence
            #ret_df.loc[rsid,"transcript"] = transcript
            #ret_df.loc[rsid,"inh"] = inh
            #ret_df.loc[rsid,"freq"] = af
            #ret_df.loc[rsid,"exon"] = exon
            #ret_df.loc[rsid,"chr_loc"] = chr_loc
            #ret_df.loc[rsid,"zyg"] = zyg
            ret_df.loc[rsid,"clin_sig_list_STR"] = clin_sig_list_STR
            ret_df.loc[rsid,"clin_sig_list"] = clin_sig_list
            ret_df.loc[rsid,"clin_sig_ct"] = clin_sig_ct
            [pubmed_dict,pmid_STR_tmp,pmid_STR] = get_publications(rsid)
            if not(pmid_STR == "-9999"):
                curr_pubmed_ctr = 0
                for key in pubmed_dict:
                    if not(curr_pubmed_ctr > pubmed_df_max_nr):
                        name = str(ret_df.loc[rsid,"gene_name"]) + " " + str(ret_df.loc[rsid,"dna_change"]) + " (" + str(ret_df.loc[rsid,"prot_change"]) + ")"
                        pubmed_list.loc[publ_ctr,"name"] = name
                        pubmed_list.loc[publ_ctr,"rsid"] = rsid
                        pubmed_list.loc[publ_ctr,"PMID"] = key
                        pubmed_list.loc[publ_ctr,"title"] = pubmed_dict[key]
                        publ_ctr = publ_ctr + 1
                        curr_pubmed_ctr = curr_pubmed_ctr + 1
            if(len(ret_df.loc[rsid,"dis_and_sig"]) < 2):
                if(pd.isnull(ret_df.loc[rsid,"dis_and_sig"])):
                    ret_df.loc[rsid,"dis_and_sig"] = []
            ret_df.loc[rsid,"rsid"] = rsid
    #for key in ret_all:
    #    print("Disease Name: " + str(key))
    #    for key2 in ret_all[key]:
    #        print("SNP ID:" + key2 + "; Gene Name: " + gene_name_dict[key2])
    #        print(ret_all[key][key2])
    #print(dict_by_snpid)
    #for snp_curr in snp_data:
        #print(snp_data[snp_curr])
        #print("SNP ID: " + str(snp_curr))
        #print(snp_data[snp_curr])
        #print(gene_name_dict[snp_curr])
        #print(dict_by_snpid[snp_curr])
        #print("Classification:\t" + snp_data[snp_curr]["clin_sig"] + "\tGene:\t" + gene_name_dict[snp_curr] + "\tDNA change:\t" + snp_data[snp_curr]["dna_change"] + "\tProtein change:\t" + snp_data[snp_curr]["prot_change"] + "\tAssociated disease(s):\t" + ";".join([i.split(";")[0] for i in dict_by_snpid[snp_curr]]))
        #for dis_and_sig in dict_by_snpid[snp_curr]:
        #    if(len(str(dis_and_sig).split(";")) > 1):
        #           print("RSID:\t" + snp_curr + "\t - \t" + str(dis_and_sig).split(";")[0] + ";\t" + str(dis_and_sig).split(";")[1] +":")
        #           if not(str(dict_by_snpid[snp_curr][dis_and_sig]) == "nan"):
        #               print("Description:\t" + str(dict_by_snpid[snp_curr][dis_and_sig]))
    comment_df = pd.DataFrame(columns=['clin_sig','name','diseasename','lab_name','comment','rsid','report_date','labrank'])
    ctr = 0
    ctr2 = 0
    comment_htmls = {}
    comment_and_disease_dict = {}
    list_of_described_genes = []
    take_list_of_comments = "true"
    for i in list(ret_df.index):
        #print(str(ret_df.loc[i,"rsid"]))
        #print(str(ret_df.loc[i,"dis_and_sig"]))
        if(isinstance(ret_df.loc[i,"dis_and_sig"], list)):
            if(len(ret_df.loc[i,"dis_and_sig"]) < 2):
                if(str(ret_df.loc[i,"rsid"]) == "nan" or pd.isnull(ret_df.loc[i,"rsid"]) or str(ret_df.loc[i,"rsid"]) == "NaN" or (str(ret_df.loc[i,"diseases"]) == "nan" or pd.isnull(ret_df.loc[i,"diseases"]) or str(ret_df.loc[i,"diseases"]) == "NaN")):
                    #print("dropped" + str(ret_df.loc[i,"rsid"]))
                    ret_df = ret_df.drop([i])
                    continue
        else:
            if(str(ret_df.loc[i,"rsid"]) == "nan" or pd.isnull(ret_df.loc[i,"rsid"]) or str(ret_df.loc[i,"rsid"]) == "NaN" or (str(ret_df.loc[i,"dis_and_sig"]) == "nan" or pd.isnull(ret_df.loc[i,"dis_and_sig"]) or str(ret_df.loc[i,"dis_and_sig"]) == "NaN") or (str(ret_df.loc[i,"diseases"]) == "nan" or pd.isnull(ret_df.loc[i,"diseases"]) or str(ret_df.loc[i,"diseases"]) == "NaN")):
                ret_df = ret_df.drop([i])
                continue
        if not(str(ret_df.loc[i,"comments_STR"]) == "nan"):
            if((len(ret_df.loc[i,"comments_STR"]) > 2) and not (str(ret_df.loc[i,"comments_STR"]) == "nan") and not (str(ret_df.loc[i,"comments_STR"])=="N/A")):
                    for k in range(0, len(list(ret_df.loc[i,"dis_and_sig"]))):
                        if not (str(list(ret_df.loc[i,"comments"])[k]) == "nan" or str(list(ret_df.loc[i,"comments"])[k]) == "N/A"):
                            name = str(ret_df.loc[i,"gene_name"]) + " " + str(ret_df.loc[i,"dna_change"]) + " (" + str(ret_df.loc[i,"prot_change"]) + ") - "
                            name = name + str(list(ret_df.loc[i,"diseases"])[k])
                            comment_df.loc[name,"clin_sig"] = str(ret_df.loc[i,"clin_sig"])
                            comment_df.loc[name,"rsid"] = str(ret_df.loc[i,"rsid"])
                            temp_comment_str = str(list(ret_df.loc[i,"comments"])[k])
                            report_date = ""
                            labname = "not known"
                            if((take_list_of_comments == "true") and ("comment:'" in temp_comment_str) and ("';date:'" in temp_comment_str)):
                                list_of_comments_tmp = temp_comment_str.split("_NEXT_COMMENT_")[:-1]
                                list_for_comment_df_comment = []
                                list_for_comment_df_report_dates = []
                                list_for_comment_df_lab_name = []
                                list_for_comment_df_labrank = []
                                list_for_comment_df_disease_name = []
                                for curr_comment_temp in list_of_comments_tmp:
                                    if not(("comment:'" in curr_comment_temp) and ("';date:'" in curr_comment_temp)):
                                        continue
                                    comment = curr_comment_temp.split("comment:'")[1].split("';date:'")[0]
                                    report_date = curr_comment_temp.split("comment:'")[1].split("';date:'")[1]
                                    labname = "not known"
                                    labrank = "1000"
                                    if("';labname:'" in curr_comment_temp.split("comment:'")[1].split("';date:'")[1]):
                                        if not("';labrank:'" in curr_comment_temp.split("comment:'")[1].split("';date:'")[1]):
                                            continue
                                        if not("';disease_name:'" in curr_comment_temp.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[1].split("';labrank:'")[1]):
                                            continue
                                        report_date = curr_comment_temp.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[0]
                                        labname = curr_comment_temp.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[1].split("';labrank:'")[0]
                                        labrank = curr_comment_temp.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[1].split("';labrank:'")[1].split("';disease_name:'")[0]
                                        diseasename = curr_comment_temp.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[1].split("';labrank:'")[1].split("';disease_name:'")[1]
                                        diseasename = diseasename[:-1]
                                    list_for_comment_df_comment.append(comment)
                                    list_for_comment_df_report_dates.append(report_date)
                                    list_for_comment_df_lab_name.append(labname)
                                    list_for_comment_df_labrank.append(labrank)
                                    list_for_comment_df_disease_name.append(diseasename)
                                comment_df.loc[name,"comment"] = list_for_comment_df_comment
                                comment_df.loc[name,"report_date"] = list_for_comment_df_report_dates
                                comment_df.loc[name,"lab_name"] = list_for_comment_df_lab_name
                                comment_df.loc[name,"labrank"] = list_for_comment_df_labrank
                                comment_df.loc[name,"diseasename"] = list_for_comment_df_disease_name
                            elif("comment:'" in temp_comment_str):
                                if("';date:'" in temp_comment_str):
                                    comment = temp_comment_str.split("comment:'")[1].split("';date:'")[0]
                                    report_date = temp_comment_str.split("comment:'")[1].split("';date:'")[1]
                                    labname = "not known"
                                    if("';labname:'" in temp_comment_str.split("comment:'")[1].split("';date:'")[1]):
                                        report_date = temp_comment_str.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[0]
                                        labname = temp_comment_str.split("comment:'")[1].split("';date:'")[1].split("';labname:'")[1]
                                        labname = labname[:-1]
                                    comment_df.loc[name,"comment"] = comment
                                    comment_df.loc[name,"report_date"] = report_date
                                    comment_df.loc[name,"lab_name"] = labname
                                    comment_df.loc[name,"labrank"] = "1000"
                                    comment_df.loc[name,"diseasename"] = ""
                                else:
                                    comment_df.loc[name,"comment"] = temp_comment_str
                                    comment_df.loc[name,"report_date"] = report_date
                                    comment_df.loc[name,"lab_name"] = labname
                                    comment_df.loc[name,"labrank"] = "1000"
                                    comment_df.loc[name,"diseasename"] = ""
                            else:
                                comment_df.loc[name,"comment"] = temp_comment_str
                                comment_df.loc[name,"report_date"] = report_date
                                comment_df.loc[name,"lab_name"] = labname
                                comment_df.loc[name,"labrank"] = "1000"
                                comment_df.loc[name,"diseasename"] = ""
                            comment_df.loc[name,"gene_name"] = str(ret_df.loc[i,"gene_name"])
        #if not(str(ret_df.loc[i,"gene_name"]) in gen_dis_dict):
        #    gen_dis_dict[str(ret_df.loc[i,"gene_name"])] = list(ret_df.loc[i,"diseases"])
        #else:
        #    gen_dis_dict[str(ret_df.loc[i,"gene_name"])].append(list(ret_df.loc[i,"diseases"]))
    
    clin_sig_all_values = ret_df["clin_sig"].tolist()
    clin_sig_all_values_STR = "".join(clin_sig_all_values)
    is_pathogenic = "false"
    if("pathogenic" in clin_sig_all_values_STR):
        is_pathogenic = "true"
    disease_list = []
    disease_dict_for_ret = {}
    for i in list(ret_df.index):
        dis_tmp_1 = ret_df.loc[i,"diseases"]
        disease_list = dis_tmp_1
        if(ret_df.loc[i,"gene_name"] in gen_dis_dict):
            dis_tmp_2 = gen_dis_dict[ret_df.loc[i,"gene_name"]]
            disease_list = dis_tmp_1 + dis_tmp_2
        [new_disease_list,syndict_2] = remove_synonyms_new(disease_list,syndict_all)
        disease_dict_for_ret[ret_df.loc[i,"gene_name"]] = new_disease_list
        syndict = syndict_2.copy()
        syndict_all.update(syndict_2)
    gene_info_ret = {}
    for gene in gene_info_dict:
        if(gene in disease_dict_for_ret):
            #gene_info_ret[gene] = [gene_info_dict[gene],disease_dict_for_ret[gene]]
            gene_info_ret[gene] = [gene_info_dict[gene],disease_dict_for_ret[gene],snp_dict[gene]]
    #for i in list(comment_df.index):
    #    gene = comment_df.loc[i,"gene_name"]
    #    if (not(gene in list_of_described_genes) and (comment_df.loc[i,"gene_name"] in gene_info_dict)):
    #        comment_htmls[ctr2] = gene_info_dict[comment_df.loc[i,"gene_name"]]
    #        if not(i in comment_and_disease_dict):
    #            comment_and_disease_dict[i] = ctr2
    #        else:
    #            comment_and_disease_dict[i].append(ctr2)
    #        ctr2 = ctr2 + 1
    #    if not(i in comment_and_disease_dict):
    #        comment_and_disease_dict[i] = ctr2
    #    else:
    #        comment_and_disease_dict[i].append(ctr2)
    #    comment_htmls[ctr2] = comment_df.loc[i,"comment"]
    return([ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic,varsome_dict])
    #return([ret_df,comment_df,pubmed_list])
    #return[ret_all,dict_by_snpid,snp_data]
    #print(df.head())
    
    

def extract_vep_data_2(dataframe,disease_filter,start,end,sample_id):
    #df = pd.read_csv(path,sep='\t')
    df = dataframe
    #print(df)
    pubmed_list = pd.DataFrame(columns=['name','PMID','title'])
    [gen_dis_dict,dis_gen_dict] = read_gene_associations(assoc_table_path)
    [gene_dis_dict_2,dis_gene_dict_2] = read_gene_associations_2(gene_disease_groups_path)
    publ_ctr = 0
    hgnc_dict = {}
    inh_dict = {}
    gene_info_dict = {}
    snp_dict = {}
    ret_df = pd.DataFrame(columns=['rsid','clin_sig','consequence','variant_class','dna_change','prot_change','exon','transcript','chr_loc','hgvs','gene_name','hgnc','freq','max_af','max_pop','vaf','gq','ad1','ad2','ad','zyg','inh','clin_sig_list_STR','clin_sig_list','clin_sig_ct','diseases','diseases_STR','disease_groups','comments','comments_STR','dis_and_sig'])
    samplename = sample_id
    #for i in range(0,len(df.index)):
    #for i in range(1400,1600):
    #for i in range(2000,3000):
    #print(dis_gene_dict_2[disease])
    #print(gene_dis_dict_2["EPHA2"])
    #print([str(i) for i in gene_dis_dict_2])
    #print([i for i in gene_dis_dict_2 if disease_filter in gene_dis_dict_2[i]])
    if(use_all_lines == "true"):
        start_line_def = 0
        stop_line_def = 20000
    if(automatic_stop_line_def == "true"):
        if(len(df.index) > 6700):
            start_line_def = 6600
            stop_line_def = 6700
        elif(len(df.index) > 100):
            start_line_def = 0
            stop_line_def = 100
    else:
        start_line_def = start
        stop_line_def = end
    clinvar_version = "ClinVar_updated_2021Aug"
    if("ClinVar_updated_2019Dec_CLNDN" in df):
        clinvar_version = "ClinVar_updated_2019Dec"
    for i in range(start_line_def,min(len(df.index),stop_line_def)):
        if ((df.loc[i,"CLIN_SIG"]=="-") and not(")" in df.loc[i,"PolyPhen"]) and (df.loc[i,clinvar_version + "_CLNDN"] == "-")):
            #print("Exited for loop at line 703: No clinical significance")
            continue
        rsid_tmp = df.loc[i,"Existing_variation"]
        temp1 = rsid_tmp.split(",")
        rsids = [i for i in temp1 if("rs" in i)]
        rsids_unique = [i for i in temp1 if(("rs" in i) and not(i in df.index))]
        if(any(rsid_curr in rsids for rsid_curr in list(df.index.values))):
            if(len(rsids_unique) > 0):
                rsid = rsids_unique[0]
            else:
                #print("Exited for loop at line 713: No RSID given")
                continue
        if(len(rsids) > 0):
            rsid = rsids[0]
            variant_class = df.loc[i,"VARIANT_CLASS"]
        else:
            #print("Exited for loop at line 719: No RSID given.")
            continue
        gene_name = df.loc[i,"SYMBOL"]
        if (not(disease_filter in dis_gene_dict_2) and (not(disease_filter == "all"))):
            #print("Disease group does not exist! Returning empty dataframes.")
            return([pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()])
        if not(disease_filter=="all"):
            if not(gene_name in gene_dis_dict_2):
                #print(str(gene_name) + "not in gene_dis_dict")
                #print("Gene name not in disease gene list")
                continue
            #else:
            #    print(gene_dis_dict_2[gene_name])
            #print(str(disease))
            elif(not(str(disease_filter) in ",".join(gene_dis_dict_2[gene_name]))):
                #print(str(disease_filter))
                #print(gene_dis_dict_2[gene_name])
                #print("Exited for loop at line 736: Gene name not in disease gene list")
                continue
        #print(i)
        hgnc_id = str(df.loc[i,"HGNC_ID"])
        if not(rsid in hgnc_dict):
            hgnc_dict[rsid] = hgnc_id
        prot_change = df.loc[i,"HGVSp"]
        if(prot_change == "-" or not("ENSP" in prot_change) or (len(prot_change.split(":")) < 2)):
            #print("Exited for loop. No Protein change given")
            continue
        prot_change = prot_change.split(":")[1]
        prot_change = prot_change.replace("%3D","=")
        consequence = df.loc[i,"Consequence"]
        consequence = consequence.split(",")[0]
        #exon = str(df.loc[i,"EXON"])
        #if(exon == "-"):
        #    exon = str(df.loc[i,"INTRON"])
        #exon = str(exon.split("/")[0])
        exon = "1"
        chr_loc = df.loc[i,"Location"]
        if not("Chr" in chr_loc or "chr" in chr_loc):
            chr_loc = "Chr"+str(chr_loc)
        af = df.loc[i,"gnomAD_AF"]
        max_af = df.loc[i,"MAX_AF"]
        if(str(max_af) == "-"):
            max_af = "N/A"
        elif(not(str(max_af).isdigit()) and str(max_af).replace('.','',1).isdigit()):
            max_af = float(max_af)
            max_af = str("{:.2%}".format(max_af))
        elif(str(max_af) == "0" or str(max_af) == "1"):
            max_af = float(max_af)
            max_af = str("{:.2%}".format(max_af))
        else:
            max_af = "N/A"
        
        max_pop = df.loc[i,"MAX_AF_POPS"]
        max_pop = max_pop.replace("gnomAD_","")
        if not(max_pop == "-"):
            max_pop_arr = max_pop.split(",")
            max_pop_arr_2 = list(dict.fromkeys(max_pop_arr))
            max_pop_arr_3 = []
            for pop in max_pop_arr_2:
                if(pop in pop_dict):
                    max_pop_arr_3.append(pop_dict[pop])
                else:
                    max_pop_arr_3.append(pop)
            max_pop = ",".join(max_pop_arr_3)
        else:
            max_pop = "N/A"
        #if("Sample_1.VAF" in df):
        vaf = df.loc[i, samplename +".VAF"]
        if not(str(vaf) == "NA"):
            if(str(max_af).replace('.','',1).isdigit()):
                vaf = float(vaf)
                vaf = str("{:.2%}".format(vaf))
        else:
            vaf = "N/A"
        #if("Sample_1.GQ" in df):
        gq = df.loc[i,samplename +".GQ"]
        #elif("GQ" in df):
        #    vaf = df.loc[i,"GQ"]
        #else:
        #    gq = 0.0
        #gq = df.loc[i,"Sample_1.GQ"]
        if(str(gq).isdigit()):
            gq = float(gq)
        ad = df.loc[i,samplename +".AD"]
        #print(str(ad))
        #print(str(ad).split(","))
        ad1 = ""
        ad2 = ""
        if(str(ad).split(",")[0].isdigit() and (len(str(ad).split(",")) > 1)):
            ad = df.loc[i,samplename +".AD"]
            ad1 = float(str(ad).split(",")[0])
            ad2 = float(str(ad).split(",")[1])
            ad = int(round(ad1 + ad2))
        else:
            ad = "NA"
            #ad = str(ad1) + "/" + str(ad2)
        zyg_tmp = df.loc[i,samplename +".GT"]
        #print(zyg_tmp)
        if(zyg_tmp.split("/")[0] == "." or (len(zyg_tmp.split("/")) < 1)):
            zyg = "---"
            #print("No zygosity given")
            continue
        elif(str(zyg_tmp.split("/")[0]) == str(zyg_tmp.split("/")[1])):
            zyg = "homozygous"
        else:
            zyg = "heterozygous"
        #print(zyg)
        if(str(af) == "-"):
            af = 999.0
        else:
            af = float(str(af))
            af = str("{:.2%}".format(af))
        #if(hgnc_id in inh_dict):
        #    inh = inh_dict[hgnc_id]
        #else:
        #    inh = get_inheritance_mode_2(hgnc_id)
        #    inh_dict[hgnc_id] = inh
        if(gene_name in inh_dict):
            inh = inh_dict[gene_name]
        else:
            inh = get_inheritance_mode_3(gene_name,str(path_prefix + "tableExport.csv"),"","false")
            inh_dict[gene_name] = inh
        #inh = get_inheritance_mode(gene_name)
        #print(inh)
        #print(prot_change)
        clin_sig_conf_tmp = df.loc[i,clinvar_version + "_CLNSIGCONF"]
        clin_sig_list_STR = ""
        clin_sig_list = []
        clin_sig_ct = []
        if not(str(clin_sig_conf_tmp) == "-"):
            sig_list = clin_sig_conf_tmp.split("%3B")
            if(len(sig_list) > 1):
                for n in sig_list:
                    clin_sig_list_STR = clin_sig_list_STR + ";" + str(n.split("(")[0])
                    clin_sig_list.append(str(n.split("(")[0]))
                    clin_sig_ct.append(str(n.split("(")[1].split(")")[0]))
            else:
                clin_sig_list_STR = clin_sig_conf_tmp.split("(")[0]
        clin_sig_2 = ""
        if(len(clin_sig_list) > 1):
            index_main_sig = clin_sig_ct.index(max(clin_sig_ct))
            if(index_main_sig < len(clin_sig_list)):
                clin_sig_2 = clin_sig_list[index_main_sig]
        clin_sig = df.loc[i,"CLIN_SIG"]
        clin_sig.replace("_"," ")
        has_clinvar = "false"
        if("," in clin_sig):
            if(len(str(df.loc[i,clinvar_version + "_CLNSIG"])) > 2):
                clin_sig = df.loc[i,clinvar_version + "_CLNSIG"]
                clin_sig = clin_sig.replace("_"," ")
        if(clin_sig == "-"):
            polyphen = df.loc[i,"PolyPhen"]
            if(")" in polyphen):
                #print("polyphen")
                #print(polyphen)
                polyphen_sig = polyphen.split("(")[0]
                polyphen_val = float(polyphen.split("(")[1].split(")")[0])
                clin_sig = polyphen_sig.replace("probably_damaging","likely pathogenic")
                clin_sig = clin_sig.replace("_"," ")
                #print(str(polyphen_val))
                #print(clin_sig)
                if not (clin_sig == "benign" or clin_sig == "likely pathogenic" or (polyphen_val > polyphen_threshold)):
                    #print("No clinical significance, polyphen value given")
                    continue
            else:
                #print("No clinical significance")
                continue
        else:
            has_clinvar = "true"
        if(clin_sig_2 == ""):
            clin_sig_2 = clin_sig
        gene_ensembl_id = df.loc[i,"Gene"]
        if("ENSG" in gene_ensembl_id and not( gene_name in gene_info_dict)):
            [gene_entrez_id,summary] = get_gene_info(gene_ensembl_id)
            gene_info_dict[gene_name] = summary
        disease = df.loc[i,clinvar_version + "_CLNDN"]
        #if(disease == "-"):
        #    continue
        if not(disease == "-"):
            disease = disease.replace("_"," ")
            diseases_curr_rs = disease.split("|")
            #print(diseases_curr_rs)
            #print(i)
            #print(diseases_curr_rs)
            diseases_curr_rs_tmp = [j for j in diseases_curr_rs if not(("not specified" in j) or ("not provided" in j) or ("not_specified" in j) or ("not_provided" in j))]
            if(len(diseases_curr_rs_tmp) > 1):
                diseases_curr_rs = remove_synonyms(diseases_curr_rs_tmp)
            else:
                diseases_curr_rs = diseases_curr_rs_tmp
        else:
            diseases_curr_rs = []
        #print(diseases_curr_rs)
        #for j in range(0,arrlen):
        #    if(("not specified" in diseases_curr_rs[j]) or ("not provided" in diseases_curr_rs[j]) or ("not_specified" in diseases_curr_rs[j]) or ("not_provided" in diseases_curr_rs[j])):
        #        diseases_curr_rs.pop(j)
        #        j = j - 1
        #        arrlen = arrlen - 1
        #print(disease)
        if(has_clinvar == "true"):
            hgvs_list = get_hgvs_list(rsid)
        else:
            hgvs_list = []
        gene_name_dict[rsid] = gene_name
        prot_change_dict[rsid] = prot_change
        snp_data_tmp = {}
        hgvsc = df.loc[i,"HGVSc"]
        if(":" in hgvsc):
            transcript = hgvsc.split(":")[0]
            dna_change = hgvsc.split(":")[1]
        else:
            transcript = "-"
            dna_change = "-"
        snp_data_tmp = {}
        if (rsid in list(df.index.values)):
            #print("RSID in dataframe index")
            continue
        if not(gene_name_dict[rsid] in snp_dict):
            snp_dict[gene_name_dict[rsid]] = []
        snp_dict[gene_name_dict[rsid]].append(rsid)
        snp_data_tmp["clin_sig"] = clin_sig
        snp_data_tmp["dna_change"] = dna_change
        snp_data_tmp["prot_change"] = prot_change
        ret_df.loc[rsid,"clin_sig"] = clin_sig
        ret_df.loc[rsid,"dna_change"] = dna_change
        ret_df.loc[rsid,"prot_change"] = prot_change
        ret_df.loc[rsid,"freq"] = af
        ret_df.loc[rsid,"max_af"] = max_af
        ret_df.loc[rsid,"max_pop"] = max_pop
        ret_df.loc[rsid,"vaf"] = vaf
        #print(gq)
        #print(rsid)
        #print(inh)
        ret_df.loc[rsid,"gq"] = gq
        ret_df.loc[rsid,"ad1"] = ad1
        ret_df.loc[rsid,"ad2"] = ad2
        ret_df.loc[rsid,"ad"] = ad
        ret_df.loc[rsid,"hgnc_id"] = hgnc_id
        ret_df.loc[rsid,"inh"] = inh
        ret_df.loc[rsid,"exon"] = exon
        ret_df.loc[rsid,"chr_loc"] = chr_loc
        ret_df.loc[rsid,"consequence"] = consequence
        ret_df.loc[rsid,"transcript"] = transcript
        ret_df.loc[rsid,"zyg"] = zyg
        ret_df.loc[rsid,"clin_sig_list_STR"] = clin_sig_list_STR
        ret_df.loc[rsid,"clin_sig_list"] = clin_sig_list
        ret_df.loc[rsid,"clin_sig_ct"] = clin_sig_ct
        ret_df.loc[rsid,"rsid"] = rsid
        ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
        [hgvs_curr,link] = get_correct_hgvs(hgvs_list,gene_name,prot_change)
        #ret = extract(link)
        if(len(hgvs_list) < 1 or (hgvs_curr == "")):
            if("benign" in clin_sig and remove_benign == "true" and not("drug response" in clin_sig)):
                #print("no clinical significance, benign mutation")
                continue
            hgvs_curr = "-"
            #print("SNP ID:" + str(rsid) + "; Gene Name: " + gene_name + "\n")
            #print("Associated disease(s):\n")
            #print(disease)
            #if(rsid in dict_by_snpid):
            #    for dis in diseases_curr_rs:
            #        dict_by_snpid[rsid][disease] = "N/A"
            #else:
            #    for dis in diseases_curr_rs:
            #        dict_by_snpid[rsid] = {}
            #        dict_by_snpid[rsid][dis] = "N/A"
            #for dis in diseases_curr_rs:
            #    if(dis in ret_all):
            #        ret_all[dis][rsid]="N/A"
            #    else:
            #        ret_all[dis] = {}
            #        ret_all[dis][rsid]="N/A"
            if(rsid in dict_by_snpid):
                for dis in diseases_curr_rs:
                    dict_by_snpid[rsid][str(dis + ";" + clin_sig)] = "N/A"
            else:
                for dis in diseases_curr_rs:
                    dict_by_snpid[rsid] = {}
                    dict_by_snpid[rsid][str(dis + ";" + clin_sig)] = "N/A"
            for dis in diseases_curr_rs:
                if(str(dis + ";" + clin_sig) in ret_all):
                    ret_all[str(dis + ";" + clin_sig)][rsid]="N/A"
                else:
                    ret_all[str(dis + ";" + clin_sig)] = {}
                    ret_all[str(dis + ";" + clin_sig)][rsid]="N/A"
            if(rsid not in dict_by_snpid):
                #print("RSID not in SNP ID dictionary")
                continue
            ret_df.loc[rsid,"dis_and_sig"] = [k for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"diseases"] = [k.split(";")[0] for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"diseases_STR"] = ";".join([k.split(";")[0] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments"] = [dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]]
            if not(gene_name in gene_dis_dict_2):
                ret_df.loc[rsid,"disease_groups"] = ""
            else:
                ret_df.loc[rsid,"disease_groups"] = ";\n".join(gene_dis_dict_2[gene_name])
            #ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
            #ret_df.loc[rsid,"chr_loc"] = chr_loc
            #ret_df.loc[rsid,"rsid"] = rsid
            #ret_df.loc[rsid,"exon"] = exon
            #ret_df.loc[rsid,"inh"] = inh
            #ret_df.loc[rsid,"freq"] = af
            #print([dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments_STR"] = ";".join([str(dict_by_snpid[rsid][k]) for k in dict_by_snpid[rsid]])
            if(len(diseases_curr_rs) > 1):
                snp_data_tmp["hgvs"] = "-"
                ret_df.loc[rsid,"hgvs"] = "-"
                snp_data[rsid]=snp_data_tmp
            else:
                #print("No diseases associated to mutation")
                continue
        else:
            #print("SNP ID:" + str(rsid) + "; Gene Name: " + gene_name + "\n")
            #print("Associated disease(s):\n")
            #[hgvs_curr,link] = get_correct_hgvs(hgvs_list,gene_name,prot_change)
            ret = extract(link)
            dis_and_sig_from_vep_keys = [(str(dis) + ";" + str(clin_sig)) for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            dis_and_sig_from_vep_values = ["N/A" for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            if("benign" in clin_sig and remove_benign=="true" and not("drug response" in clin_sig)):
                dis_and_sig_from_vep_keys = []
                dis_and_sig_from_vep_values = []
            dis_and_sig_curr = dict(zip(dis_and_sig_from_vep_keys,dis_and_sig_from_vep_values))
            #dis_and_sig_curr_rs = [(str(dis) + ";" + str(clin_sig)) for dis in diseases_curr_rs if not(("not specified" in dis) or ("not provided" in dis) or ("not_specified" in dis) or ("not_provided" in dis))]
            #ret2 = remove_synonyms_2[ret]
            #print("lists: ")
            #print(dis_and_sig_from_vep_keys)
            #print(list(ret))
            for key in list(ret):
                if(("not specified" in key) or ("not provided" in key) or ("not_specified" in key) or ("not_provided" in key)):
                    del ret[key]
                elif("benign" in key and remove_benign=="true" and not("drug response" in clin_sig)):
                    del ret[key]
                else:
                    dis_and_sig_curr[key] = ret[key]
                    #dis_and_sig_curr_rs.append(key)
            #print(i)
            dis_and_sig_curr_rs = remove_synonyms_2(dis_and_sig_curr)
            #print(dis_and_sig_curr_rs)
            #diseases_curr_rs_tmp = [j for j in diseases_curr_rs if not(("not specified" in j) or ("not provided" in j) or ("not_specified" in j) or ("not_provided" in j))]
            #print("d.split(;)[0]")
            #print([d.split(";")[0] for d in dis_and_sig_curr_rs])
            #dis_no_syn = remove_synonyms([d.split(";")[0] for d in dis_and_sig_curr_rs])
            #print(dis_no_syn)
            #dis_and_sig_curr_rs = [dis for dis in dis_and_sig_curr_rs if (dis.split(";")[0] in dis_no_syn)]
            if(len(ret) < 1):
                continue
            for dis in dis_and_sig_curr_rs:
                if(rsid in dict_by_snpid):
                    dict_by_snpid[rsid][dis] = dis_and_sig_curr_rs[dis]
                    #dict_by_snpid[rsid][dis] = ret[dis]
                else:
                    dict_by_snpid[rsid] = {}
                    dict_by_snpid[rsid][dis] = dis_and_sig_curr_rs[dis]
                    #dict_by_snpid[rsid][dis] = ret[dis]
                if(dis in ret_all):
                    ret_all[dis][rsid]= dis_and_sig_curr_rs[dis]
                    #ret_all[dis][rsid]=ret[dis]
                else:
                    ret_all[dis] = {}
                    ret_all[dis][rsid]= dis_and_sig_curr_rs[dis]
                    #ret_all[dis][rsid]=ret[dis]
                if(dis not in dis_and_sig_from_vep_keys):
                    if(rsid not in dict_by_snpid):
                        #print("RSID not in SNP dictionary")
                        continue
            #print(dict_by_snpid[rsid])
            #snp_data_tmp[rsid]=dict_by_snpid[rsid]
            snp_data_tmp["hgvs"] = hgvs_curr
            ret_df.loc[rsid,"hgvs"] = hgvs_curr
            snp_data_tmp["gene_name"] = gene_name
            ret_df.loc[rsid,"gene_name"] = gene_name
            #snp_data_tmp["diseases"] = ";".join([key for key in ret])
            snp_data[rsid]=snp_data_tmp
            #print(ret_df)
            ret_df.loc[rsid,"diseases"] = [k.split(";")[0] for k in dict_by_snpid[rsid]]
            #disease_with_inh = ""
            #for disease_name_curr in dict_by_snpid[rsid]:
            #    inh = get_inheritance_mode_3(gene_name,"VEP_test/tableExport.csv",disease_name_curr,"false")
            ret_df.loc[rsid,"diseases_STR"] = ";".join([k.split(";")[0] for k in dict_by_snpid[rsid]])
            if not(gene_name in gene_dis_dict_2):
                ret_df.loc[rsid,"disease_groups"] = ""
            else:
                ret_df.loc[rsid,"disease_groups"] = ";\n".join(gene_dis_dict_2[gene_name])
            ret_df.loc[rsid,"comments"] = [dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]]
            #print([dict_by_snpid[rsid][k] for k in dict_by_snpid[rsid]])
            ret_df.loc[rsid,"comments_STR"] = ";".join([str(dict_by_snpid[rsid][k]) for k in dict_by_snpid[rsid]])
            #ret_df.loc[rsid,"gene_name"] = gene_name_dict[rsid]
            ret_df.loc[rsid,"dis_and_sig"] = [k for k in dict_by_snpid[rsid]]
            ret_df.loc[rsid,"consequence"] = consequence
            #ret_df.loc[rsid,"transcript"] = transcript
            #ret_df.loc[rsid,"inh"] = inh
            #ret_df.loc[rsid,"freq"] = af
            #ret_df.loc[rsid,"exon"] = exon
            #ret_df.loc[rsid,"chr_loc"] = chr_loc
            #ret_df.loc[rsid,"zyg"] = zyg
            ret_df.loc[rsid,"clin_sig_list_STR"] = clin_sig_list_STR
            ret_df.loc[rsid,"clin_sig_list"] = clin_sig_list
            ret_df.loc[rsid,"clin_sig_ct"] = clin_sig_ct
            [pubmed_dict,pmid_STR_tmp,pmid_STR] = get_publications(rsid)
            if not(pmid_STR == "-9999"):
                for key in pubmed_dict:
                    name = str(ret_df.loc[rsid,"gene_name"]) + " " + str(ret_df.loc[rsid,"dna_change"]) + " (" + str(ret_df.loc[rsid,"prot_change"]) + ")"
                    pubmed_list.loc[publ_ctr,"name"] = name
                    pubmed_list.loc[publ_ctr,"rsid"] = rsid
                    pubmed_list.loc[publ_ctr,"PMID"] = key
                    pubmed_list.loc[publ_ctr,"title"] = pubmed_dict[key]
                    publ_ctr = publ_ctr + 1
            if(len(ret_df.loc[rsid,"dis_and_sig"]) < 2):
                if(pd.isnull(ret_df.loc[rsid,"dis_and_sig"])):
                    ret_df.loc[rsid,"dis_and_sig"] = []
            ret_df.loc[rsid,"rsid"] = rsid
    #for key in ret_all:
    #    print("Disease Name: " + str(key))
    #    for key2 in ret_all[key]:
    #        print("SNP ID:" + key2 + "; Gene Name: " + gene_name_dict[key2])
    #        print(ret_all[key][key2])
    #print(dict_by_snpid)
    #for snp_curr in snp_data:
        #print(snp_data[snp_curr])
        #print("SNP ID: " + str(snp_curr))
        #print(snp_data[snp_curr])
        #print(gene_name_dict[snp_curr])
        #print(dict_by_snpid[snp_curr])
        #print("Classification:\t" + snp_data[snp_curr]["clin_sig"] + "\tGene:\t" + gene_name_dict[snp_curr] + "\tDNA change:\t" + snp_data[snp_curr]["dna_change"] + "\tProtein change:\t" + snp_data[snp_curr]["prot_change"] + "\tAssociated disease(s):\t" + ";".join([i.split(";")[0] for i in dict_by_snpid[snp_curr]]))
        #for dis_and_sig in dict_by_snpid[snp_curr]:
        #    if(len(str(dis_and_sig).split(";")) > 1):
        #           print("RSID:\t" + snp_curr + "\t - \t" + str(dis_and_sig).split(";")[0] + ";\t" + str(dis_and_sig).split(";")[1] +":")
        #           if not(str(dict_by_snpid[snp_curr][dis_and_sig]) == "nan"):
        #               print("Description:\t" + str(dict_by_snpid[snp_curr][dis_and_sig]))
    comment_df = pd.DataFrame(columns=['clin_sig','name','comment','rsid'])
    ctr = 0
    ctr2 = 0
    comment_htmls = {}
    comment_and_disease_dict = {}
    list_of_described_genes = []
    for i in list(ret_df.index):
        #print(str(ret_df.loc[i,"rsid"]))
        #print(str(ret_df.loc[i,"dis_and_sig"]))
        if(isinstance(ret_df.loc[i,"dis_and_sig"], list)):
            if(len(ret_df.loc[i,"dis_and_sig"]) < 2):
                if(str(ret_df.loc[i,"rsid"]) == "nan" or pd.isnull(ret_df.loc[i,"rsid"]) or str(ret_df.loc[i,"rsid"]) == "NaN" or (str(ret_df.loc[i,"diseases"]) == "nan" or pd.isnull(ret_df.loc[i,"diseases"]) or str(ret_df.loc[i,"diseases"]) == "NaN")):
                    #print("dropped" + str(ret_df.loc[i,"rsid"]))
                    ret_df = ret_df.drop([i])
                    #print("Exited for loop at line 1161: RSID is NaN or null")
                    continue
        else:
            if(str(ret_df.loc[i,"rsid"]) == "nan" or pd.isnull(ret_df.loc[i,"rsid"]) or str(ret_df.loc[i,"rsid"]) == "NaN" or (str(ret_df.loc[i,"dis_and_sig"]) == "nan" or pd.isnull(ret_df.loc[i,"dis_and_sig"]) or str(ret_df.loc[i,"dis_and_sig"]) == "NaN") or (str(ret_df.loc[i,"diseases"]) == "nan" or pd.isnull(ret_df.loc[i,"diseases"]) or str(ret_df.loc[i,"diseases"]) == "NaN")):
                ret_df = ret_df.drop([i])
                #print("RSID is null or NaN or diseases are null or NaN")
                continue
        if not(str(ret_df.loc[i,"comments_STR"]) == "nan"):
            if((len(ret_df.loc[i,"comments_STR"]) > 2) and not (str(ret_df.loc[i,"comments_STR"]) == "nan") and not (str(ret_df.loc[i,"comments_STR"])=="N/A")):
                    for k in range(0, len(list(ret_df.loc[i,"dis_and_sig"]))):
                        if not (str(list(ret_df.loc[i,"comments"])[k]) == "nan" or str(list(ret_df.loc[i,"comments"])[k]) == "N/A"):
                            name = str(ret_df.loc[i,"gene_name"]) + " " + str(ret_df.loc[i,"dna_change"]) + " (" + str(ret_df.loc[i,"prot_change"]) + ") - "
                            name = name + str(list(ret_df.loc[i,"diseases"])[k])
                            comment_df.loc[name,"clin_sig"] = str(ret_df.loc[i,"clin_sig"])
                            comment_df.loc[name,"rsid"] = str(ret_df.loc[i,"rsid"])
                            comment_df.loc[name,"comment"] = str(list(ret_df.loc[i,"comments"])[k])
                            comment_df.loc[name,"gene_name"] = str(ret_df.loc[i,"gene_name"])
        #if not(str(ret_df.loc[i,"gene_name"]) in gen_dis_dict):
        #    gen_dis_dict[str(ret_df.loc[i,"gene_name"])] = list(ret_df.loc[i,"diseases"])
        #else:
        #    gen_dis_dict[str(ret_df.loc[i,"gene_name"])].append(list(ret_df.loc[i,"diseases"]))
    clin_sig_all_values = ret_df["clin_sig"].tolist()
    clin_sig_all_values_STR = "".join(clin_sig_all_values)
    is_pathogenic = "false"
    if("pathogenic" in clin_sig_all_values_STR):
        is_pathogenic = "true"
    disease_list = []
    disease_dict_for_ret = {}
    for i in list(ret_df.index):
        dis_tmp_1 = ret_df.loc[i,"diseases"]
        disease_list = dis_tmp_1
        if(ret_df.loc[i,"gene_name"] in gen_dis_dict):
            dis_tmp_2 = gen_dis_dict[ret_df.loc[i,"gene_name"]]
            disease_list = dis_tmp_1 + dis_tmp_2
        disease_dict_for_ret[ret_df.loc[i,"gene_name"]] = remove_synonyms(disease_list)
    gene_info_ret = {}
    for gene in gene_info_dict:
        if(gene in disease_dict_for_ret):
            #gene_info_ret[gene] = [gene_info_dict[gene],disease_dict_for_ret[gene]]
            gene_info_ret[gene] = [gene_info_dict[gene],disease_dict_for_ret[gene],snp_dict[gene]]
    #for i in list(comment_df.index):
    #    gene = comment_df.loc[i,"gene_name"]
    #    if (not(gene in list_of_described_genes) and (comment_df.loc[i,"gene_name"] in gene_info_dict)):
    #        comment_htmls[ctr2] = gene_info_dict[comment_df.loc[i,"gene_name"]]
    #        if not(i in comment_and_disease_dict):
    #            comment_and_disease_dict[i] = ctr2
    #        else:
    #            comment_and_disease_dict[i].append(ctr2)
    #        ctr2 = ctr2 + 1
    #    if not(i in comment_and_disease_dict):
    #        comment_and_disease_dict[i] = ctr2
    #    else:
    #        comment_and_disease_dict[i].append(ctr2)
    #    comment_htmls[ctr2] = comment_df.loc[i,"comment"]
    return([ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic])
    #return([ret_df,comment_df,pubmed_list])
    #return[ret_all,dict_by_snpid,snp_data]
    #print(df.head())

def get_gene_info(ensg):
    gene_info = mg.getgene(ensg, fields=['name'])
    if gene_info is None:
        return(["",""])
    if("summary" in gene_info):
        summary = gene_info["summary"]
    else:
        summary = ""
    if("_id" in gene_info):
        gene_id = gene_info["_id"]
    else:
        gene_id = ""
    return([gene_id,summary])

def extract(url):
    #url = sys.argv[1]
    #print(url)
    response = urllib.request.urlopen(url)
    webContent = response.read()
    dfs = pd.read_html(webContent)
    if(len(dfs) < 3):
        return ({})
    table_as_df = dfs[2]
    submitters = []
    significances = []
    significance_categories = []
    conditions = []
    comments = []
    dates = []
    give_list_of_comments = "true"
    infile = open(file_with_labranks)
    labrank_string=infile.read()
    infile.close()
    list_with_labranks = labrank_string.split("\n")
    ranking = {}
    ctr = 1
    for lab_curr in list_with_labranks:
        ranking[lab_curr] = ctr
        ctr = ctr + 1
    #ranking["ClinGen RASopathy Variant Curation Expert Panel"] = 10
    #ranking["Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine"] = 9
    #ranking["GeneDx"] = 8
    #ranking["Invitae"] = 1
    #ranking["Mendelics"] = 11
    list_of_diseases = []
    name_and_significances = []
    cond_ranking = []
    for i in range(0,len(table_as_df.index)):
        date = table_as_df.iloc[i,5]
        list_of_conditions = table_as_df.iloc[i,4].split(";")
        for j in range(0,len(list_of_conditions)):
            submitters.append(table_as_df.iloc[i,0])
            if not(table_as_df.iloc[i,0] in ranking):
                cond_ranking.append(99)
            else:
                cond_ranking.append(ranking[table_as_df.iloc[i,0]])
            significances.append(table_as_df.iloc[i,3])
            conditions.append(list_of_conditions[j])
            comments.append(table_as_df.iloc[i,8])
            sig_tmp = "neu"
            if((table_as_df.iloc[i,3] == "pathogenic") or (table_as_df.iloc[i,3] == "likely pathogenic")):
                sig_tmp = "pat"
            if((table_as_df.iloc[i,3] == "benign") or (table_as_df.iloc[i,3] == "likely benign")):
                sig_tmp = "ben"
            name_and_significances.append(str(table_as_df.iloc[i,0]) + ";" + sig_tmp)
            significance_categories.append(sig_tmp)
            if not(list_of_conditions[j] in list_of_diseases):
                list_of_diseases.append(list_of_conditions[j])
            dates.append(date)
    significance_ranking = {}
    significance_ranking["pathogenic"] = 2
    significance_ranking["likely pathogenic"] = 1
    significance_ranking["benign"] = 2
    significance_ranking["likely benign"] = 1
    #cond_ranking = [ranking[i] for i in submitters]
    for i in range(0,len(comments)):
        if not((str(significances[i]) == "pathogenic") or (str(significances[i]) == "likely pathogenic")):
            cond_ranking[i] = 999
        if not(str(significances[i]) in ranking):
            cond_ranking[i] = 999
        if((not comments[i]) or (len(str(comments[i]))<4) or (str(comments[i]) == "")):
            cond_ranking[i] = 999
    ret = {}
    dict_of_comments = {}
    if(sort_by_clin_sig == "false"):
        for dis in list_of_diseases:
            for sig in ["pat","ben"]:
                #ranking_temp = [cond_ranking[i] for i in range(0,len(cond_ranking)) if (str(conditions[i]) == str(dis) and str(significance_categories[i]) == sig)]
                #this_disease = [i for i in range(0,len(cond_ranking)) if (str(conditions[i]) == str(dis))]
                #this_disease = [i for i in range(0,len(cond_ranking)) if (str(conditions[i]) == str(dis) and str(significance_categories[i]) == sig)]
                ranking_temp = [cond_ranking[i] for i in range(0,len(cond_ranking)) if (str(dis) in str(conditions[i]) and str(significance_categories[i]) == sig)]
                this_disease = [i for i in range(0,len(cond_ranking)) if (str(dis) in str(conditions[i]) and str(significance_categories[i]) == sig)]
                if(len(ranking_temp) > 0):
                    if(give_list_of_comments == "true"):
                        sig_new = ""
                        if(sig == "pat"):
                            sig_new = "pathogenic"
                        else:
                            sig_new = "benign"
                        ret[dis + ";" + sig_new] = ""
                        for i in range(0,len(ranking_temp)):
                            if not(int(ranking_temp[i]) in dict_of_comments):
                                dict_of_comments[int(ranking_temp[i])] = []
                            curr_comment = {}
                            curr_text = comments[i]
                            dis_and_sig = dis + ";" + significances[i]
                            #curr_comment[str(dis_and_sig)] = "comment:'" + str(curr_text) + "';date:'" + str(dates[i]) + "';labname:'" + str(submitters[i]) + "';labrank:'" + str(ranking_temp[i]) + "'"
                            labrank = 999
                            if(submitters[i] in ranking):
                                labrank = ranking[submitters[i]]
                            curr_comment[str(dis + ";" + sig_new)] = "comment:'" + str(curr_text) + "';date:'" + str(dates[i]) + "';labname:'" + str(submitters[i]) + "';labrank:'" + str(labrank) + "';disease_name:'" + str(conditions[i]) + "'"
                            ret[dis + ";" + sig_new] = ret[dis + ";" + sig_new] + curr_comment[str(dis + ";" + sig_new)] + "_NEXT_COMMENT_"
                            if(int(ranking_temp[i]) in dict_of_comments):
                                dict_of_comments[int(ranking_temp[i])].append(curr_comment)
                            else:
                                dict_of_comments[int(ranking_temp[i])].append(curr_comment)
                    else:
                        nbr_best_entry_this_disease = np.argmin(ranking_temp)
                        best_text = comments[nbr_best_entry_this_disease]
                        dis_and_sig = dis + ";" + significances[nbr_best_entry_this_disease]
                        ret[dis_and_sig] = "comment:'" + best_text + "';date:'" + dates[nbr_best_entry_this_disease] + "';labname:'" + submitters[nbr_best_entry_this_disease] + "'"
    else:
        for dis in list_of_diseases:
            for sig in ["pat","ben"]:
                #ranking_temp = [cond_ranking[i] for i in range(0,len(cond_ranking)) if (str(conditions[i]) == str(dis) and str(significance_categories[i]) == sig)]
                #this_disease = [i for i in range(0,len(cond_ranking)) if (str(conditions[i]) == str(dis) and str(significance_categories[i]) == sig)]
                ranking_temp = [cond_ranking[i] for i in range(0,len(cond_ranking)) if (str(dis) in str(conditions[i]) and str(significance_categories[i]) == sig)]
                this_disease = [i for i in range(0,len(cond_ranking)) if (str(dis) in str(conditions[i]) and str(significance_categories[i]) == sig)]
                significance_ranking_tmp = [significance_ranking[significances[this_disease[i]]] for i in range(0,len(this_disease)) if    (str(significance_categories[this_disease[i]]) == sig)]
                if(len(significance_ranking_tmp) > 0):
                    best_significance_tmp = min(significance_ranking_tmp)
                    ranking_temp_2 = [ranking_temp[i] for i in range(0,len(ranking_temp)) if (significance_ranking_tmp[i] == best_significance_tmp)]
                    this_disease_2 = [i for i in range(0,len(ranking_temp)) if (significance_ranking_tmp[i] == best_significance_tmp)]
                    if(len(ranking_temp_2) > 0):
                        nbr_best_entry_this_disease = this_disease[np.argmin(ranking_temp_2)]
                        best_text = comments[nbr_best_entry_this_disease]
                        dis_and_sig = dis + ";" + significances[nbr_best_entry_this_disease]
                        ret[dis_and_sig] = "comment:'" + best_text + "';date:'" + dates[nbr_best_entry_this_disease] + "';labname:'" + submitters[nbr_best_entry_this_disease] + "'"
                        #ret[dis] = best_text
    if(give_list_of_comments == "true"):
        return ret
        #return dict_of_comments
    return ret


def generate_report_and_return_pdf(path,disease_list,polygen_list):
    for disease_curr in disease_list:
        [ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic] = extract_vep_data(path,disease_list)
        [ret_df_copy,comment_df_copy,pubmed_list_copy] = [ret_df.copy(deep=True),comment_df.copy(deep=True),pubmed_list.copy(deep=True)]
        [table_1_html,comment_table_html,table_2_html,pubmed_html,comments_style,pubmed_style] = make_table(ret_df,comment_df,pubmed_list,gene_info_ret,disease_curr)
        filename1 = "output" + str(disease_curr.replace(" ","").replace("/","")) + ".html" 
        htmlstr=generate_disease_report("dis_calc/templates/dis_calc/disease_report_to_fill.html",disease_curr,disease_curr,filename1,table_1_html,table_2_html,comment_table_html,pubmed_html,comments_style,pubmed_style,"Max Mustermann","1.1.1900",is_pathogenic)
        filename2 = "report" + str(disease_curr.replace(" ","").replace("/","")) + ".pdf" 
        make_pdf_from_str(htmlstr,filename2)

#link = get_correct_hgvs_link(get_hgvs_list("rs28897727"),"BRCA2","p.Asp1420Tyr")
#print(link)
#ret = extract(link)
#for key in ret:
#    print(str(key) + ":\t" + str(ret[key]) + "\n")
#inh = get_inheritance_mode("AGRN")
#print(inh)
#print(get_publications("rs17160775"))
#print(get_synonyms("Idiopathic generalized epilepsy",""))
#[dict_by_disease_output,dict_by_snpid_output,snp_data_output] = extract_vep_data("/media/quirin/INTENSO/VEP_test/first_20000_lines.txt")
#[ret_df,comment_df] = extract_vep_data("/media/quirin/INTENSO/VEP_test/new_20000.txt")
#generate_json("Likely pathogenic sequence variant(s) detected .","yellow",ret_df)
#for key in ret:
#    print(str(key) + ":\t" + str(ret[key]) + "\n")
#print(ret)
#print(get_inheritance_mode_2("14825"))
#print(get_link("23336",""))
#make_pdf("VEP_test/report2.html","output_opth_report_3.pdf")

#make_pdf("VEP_test/output_ophtalmological_test.html","output_opth_report_2.pdf")
#read_gene_associations(assoc_table_path)
#print(get_gene_info("EPHA2"))
#[gene_dis_dict,dis_gene_dict] = read_gene_associations_2("/home/quirin/Downloads/gene_disease_groups.csv")
#print(dis_gene_dict)
#[ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic] = extract_vep_data("VEP_test/first_20000_lines.txt","Ophthalmological conditions")
#print(is_pathogenic)
#[ret_df_copy,comment_df_copy,pubmed_list_copy] = [ret_df.copy(deep=True),comment_df.copy(deep=True),pubmed_list.copy(deep=True)]
#[table_1_html,comment_table_html,table_2_html,pubmed_html,comments_style,pubmed_style] = make_table(ret_df,comment_df,pubmed_list,gene_info_ret,"Ophthalmological conditions")
#print(table_1_html)
#htmlstr=generate_disease_report("VEP_test/disease_report_to_fill.html","Ophthalmological conditions","Ophthalmological conditions","VEP_test/output_ophtalmological_test.html",table_1_html,table_2_html,comment_table_html,pubmed_html,comments_style,pubmed_style,"Max Mustermann","1.1.1900",is_pathogenic)
#generate_disease_report("VEP_test/disease_report_to_fill.html","Ophthalmological conditions","Ophthalmological conditions","VEP_test/output_ophtalmological_test.html",table_1_html,table_2_html,comment_table_html,pubmed_html,comments_style,pubmed_style,"Max Mustermann","1.1.1900",is_pathogenic)
#htmlstr=make_html("VEP_test/table_to_fill.html","VEP_test/output_test.html","Likely pathogenic sequence variant(s) in gene related to reported phenotype detected .",table_1_html,table_2_html,comment_table_html,pubmed_html,comments_style,pubmed_style,"Max Mustermann","1.1.1900",is_pathogenic)
#make_pdf_from_str(htmlstr,"report_ophtalmological.pdf")

#generate_json("Likely pathogenic sequence variant(s) detected .","yellow",ret_df_copy,comment_df_copy,pubmed_list_copy)

def run_all(vcf_df,diseases,outfile,sample_id,customer_id):
    report_paths = []
    if(diseases[0] == "all"):
        diseases_to_filter = ["all"]
    diseases_to_filter=[i.replace("_"," ") for i in diseases]
    #startpos1 = 6654
    #stoppos1 = 6656
    #startpos1 = 6650
    #stoppos1 = 6670
    startpos1 = 0
    stoppos1 = 6700
    cancer_groups_for_filter = ["BREAST","GASTROINTESTINAL","LUNG","PROSTATE","GENITOURINARY","CNS","FEMALE_REPRODUCTIVE","MALE_REPRODUCTIVE","SKIN","ENDOCRINE","HEMATOLOGIC","MUSCULOSKELETAL"]
    # iterate over selected disease types
    for dis in diseases_to_filter:
        outfile_curr = outfile + ".html"
        patient_id = outfile.replace("/scripts/dis_report_","")
        ## extract data from VEP output
        #[ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic] = extract_vep_data(vcf_df,dis,startpos1,stoppos1,sample_id)
        [ret_df,comment_df,pubmed_list,gene_info_ret,is_pathogenic,varsome_dict] = extract_vep_data(vcf_df,dis,startpos1,stoppos1,sample_id)
        ret_df_old = ret_df.copy()
        comment_df_old = comment_df.copy()
        pubmed_list_old = pubmed_list.copy()
        [results_json,comments_json,freqs_json,pubmed_json] = generate_json_2("none","yellow",ret_df,comment_df,pubmed_list,"false",[],customer_id,varsome_dict)
        if(dis == "all"):
            [results_json_2,comments_json_2,freqs_json_2,pubmed_json_2] = generate_json_2("none","yellow",ret_df_old,comment_df_old,pubmed_list_old,"true",cancer_groups_for_filter,customer_id,varsome_dict)
        #print("json results:")
        #print(results_json)
        #write_json_to_one_file([results_json,comments_json,freqs_json,pubmed_json],"/scripts/dis_genes_" + patient_id + "_" + dis + ".json")
        write_json_to_one_file([results_json],str(path_prefix + "dis_genes_" + patient_id + "_" + dis + ".json"))
        if(dis == "all"):
            write_json_to_one_file([results_json_2],str(path_prefix + "dis_genes_" + patient_id + "_" + dis + "_by_cancer_group.json"))
        #[ret_df_copy,comment_df_copy,pubmed_list_copy] = [ret_df.copy(deep=True),comment_df.copy(deep=True),pubmed_list.copy(deep=True)]
        # write tables with output data
        #[table_1_html,comment_table_html,table_2_html,pubmed_html,comments_style,pubmed_style] = make_table(ret_df,comment_df,pubmed_list,gene_info_ret,dis)
        #htmlstr=generate_disease_report("dis_calc/disease_report_to_fill.html",dis,dis,(outfile+".html"),table_1_html,table_2_html,comment_table_html,pubmed_html,comments_style,pubmed_style,"Max Mustermann","1.1.1900",is_pathogenic)
        # write disease report
        #htmlstr=generate_disease_report("dis_calc/disease_report_to_fill.html",dis,dis,outfile_curr,table_1_html,table_2_html,comment_table_html,pubmed_html,comments_style,pubmed_style,"Max Mustermann","1.1.1900",is_pathogenic)
        #make_pdf_from_str(htmlstr,(outfile + "_" + dis + ".pdf"))
        report_paths.append(str(path_prefix + "dis_genes_" + patient_id + "_" + dis + ".json"))
        if(dis == "all"):
            report_paths.append(str(path_prefix + "dis_genes_" + patient_id + "_" + dis + "_by_cancer_group.json"))
    return(report_paths)


#read_varsome()
print(read_varsome("rs28928907"))
#dis_and_sig = ['essential thrombocytemia;Pathogenic', 'Congenital amegakaryocytic thrombocytopenia;pathogenic', 'Thrombocytopenia;pathogenic', ' essential thrombocytemia;pathogenic']
#ret = remove_str_from_dis_and_sig("dis_str_to_remove.txt",dis_and_sig)
#print(ret)
sys.exit("")

#ret = extract("https://clinvarminer.genetics.utah.edu/submissions-by-variant/NM_006205.3%28PDE6H%29%3Ac.35C%3EG%20%28p.Ser12Ter%29")
#ret = extract("https://clinvarminer.genetics.utah.edu/submissions-by-variant/NM_005373.2%28MPL%29%3Ac.305G%3EC%20%28p.Arg102Pro%29")
#print(ret)
#sys.exit("")


if(run_as_script == "true"):
    if(len(sys.argv) < 5):
        #print("Not enough parameters")
        quit()
    patient_vcf = sys.argv[1]
    sample_id = "Sample_1"
    if(len(patient_vcf.split("_")) > 1):
        path_new = patient_vcf.split("/")[len(patient_vcf.split("/")) -1]
        if(path_new == "file_for_disease_genes.txt"):
            sample_id = "Sample_1"
        elif(path_new == "first_20000_lines.txt"):
            sample_id = "Sample_1"
        else:
            sample_id = path_new.split("_")[0]
    disease_group_file = sys.argv[2]
    outfile = sys.argv[3]
    customer_id = sys.argv[4]
    if(disease_group_file == "all"):
        diseases = ["all"]
    else:
        disgr_file = open(disease_group_file)
        disgr_str=disgr_file.read()
        disgr_file.close()
        lines = disgr_str.split("\n")
        diseases = []
        for line in lines:
            lineSplit = line.split("\t")
            disgr_name = lineSplit[0].replace(" ","_")
            diseases.append(disgr_name)
        diseases = ["all"]
    vcf_df = pd.read_csv(patient_vcf,sep='\t')
    report_paths = run_all(vcf_df,diseases,outfile,sample_id,customer_id)
    print("_SEPERATOR_" + ("_SEPERATOR_".join(report_paths)))
    
    
    

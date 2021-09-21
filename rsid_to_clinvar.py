import sys
import urllib.request
import urllib.error
import urllib.parse
import pandas as pd
import numpy as np
import requests
#from extract_table import extract
import os
from urllib.parse import urlparse
from urllib.request import Request
#import grequests
import mygene


snpid="rs28897727"
allele = "T"
gene_name = "BRCA2"
prot_change="p.Asp1420Tyr"


def get_ticket_OLD(tgt):
    tgt = "TGT-724267-NwRf7rf4Tjr2dO90eOFpllxu4G5kQIXvWZlvyyTHtoyd1PQ5K4-cas"
    cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/tickets/{" + tgt + "} -H 'content-type: application/x-www-form-urlencoded' -d service=http%3A%2F%2Fumlsks.nlm.nih.gov"
    stream = os.popen(cmd)
    output = stream.read()
    #print(output)
    return(output)


def get_tgt(api_key):
    cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/api-key -H 'content-type: application/x-www-form-urlencoded' -d apikey=b737cbc4-4006-44a1-b7c8-70e02102a7bd"
    stream = os.popen(cmd)
    output = stream.read()
    #print(output)
    if ("form action=\"" in output):
        strSplit1 = output.split("form action=\"")[1]
        #print(strSplit1)
        if("\"" in strSplit1):
            strSplit2 = strSplit1.split("\"")
            #print(strSplit2)
            if("api-key/" in strSplit2[0]):
                ret = strSplit2[0].split("api-key/")[1]
                return(ret)
    return("none")

def get_ticket(tgt):
    #tgt = "TGT-724267-NwRf7rf4Tjr2dO90eOFpllxu4G5kQIXvWZlvyyTHtoyd1PQ5K4-cas"
    tgt = "TGT-233081-Oo3Bxdw1x7LuCRnmJHTXuiPQfJae0esgaMnSbGU6ccN1geRfbb-cas"
    
    #cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/tickets/{" + tgt + "} -H 'content-type: application/x-www-form-urlencoded' -d service=http%3A%2F%2Fumlsks.nlm.nih.gov"
    
    cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/api-key/" + tgt + " -H 'content-type: application/x-www-form-urlencoded' -d service=http%3A%2F%2Fumlsks.nlm.nih.gov"
    #cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/api-key/TGT-243087-XyxjqTP91rsLiiMetMLQulxugY6kDaSVCQSBym7F1vQnFtmSsF-cas -H 'content-type: application/x-www-form-urlencoded' -d service=http%3A%2F%2Fumlsks.nlm.nih.gov"
    stream = os.popen(cmd)
    output = stream.read()
    if(len(output) < 1):
        tgt_2 = get_tgt("")
        if not(tgt_2 == "none"):
            tgt = tgt_2
            cmd = "curl -X POST https://utslogin.nlm.nih.gov/cas/v1/api-key/" + tgt + " -H 'content-type: application/x-www-form-urlencoded' -d service=http%3A%2F%2Fumlsks.nlm.nih.gov"
            stream = os.popen(cmd)
            output = stream.read()
    #print(output)
    #print(output)
    return(output)

def is_synonym(disease_1,disease_2,tgt):
    syn = get_synonyms(disease_1,"TGT-724267-NwRf7rf4Tjr2dO90eOFpllxu4G5kQIXvWZlvyyTHtoyd1PQ5K4-cas")
    if(disease_2 in syn):
        return("true")
    else:
        return("false")

    
def get_synonyms(disease_STR,tgt):
    if("not specified" in disease_STR):
        return([])
    auth = get_ticket("TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas")
    #auth = "ST-2227617-atRsjZqEFS6te6eP1Vme-cas"
    disease = urllib.parse.quote_plus(disease_STR)
    #disease = tmp.geturl()
    #print(disease)
    url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=" + disease + "&ticket=" + auth
    request = requests.get(url)
    if not(request.status_code == 200):
        return([])
    try:
        response = urllib.request.urlopen(Request(url, headers={'User-Agent': 'Mozilla/5.0'}))
    except:
        return([])
    #req = Request('http://www.cmegroup.com/trading/products/#sortField=oi&sortAsc=false&venues=3&page=1&cleared=1&group=1', headers={'User-Agent': 'Mozilla/5.0'})
    #webpage = urlopen(req).read()
    if("403" in str(response.getcode())):
        return([])
    #table_file = open("VEP_test/test_syn.json")
    #webContent=table_file.read()
    if not("\"results\":[" in webContent):
        return([])
    result_tmp = webContent.split("\"results\":[")[1]
    result = result_tmp.split("]")[0]
    result_list = result.split("}")
    ret = []
    for i in result_list:
        if not("\"name\":\"" in i):
            continue
        resultSplit = i.split("\"name\":\"")[1]
        resultStr = resultSplit.split("\"")[0]
        ret.append(resultStr.lower())
    return(ret)


def get_synonyms_new(disease_STR,tgt,syndict):
    syndict_2 = syndict.copy()
    #print("disease name")
    #print(disease_STR)
    if("not specified" in disease_STR):
        syndict = syndict_2.copy()
        return([],syndict)
    if("not provided" in disease_STR):
        syndict = syndict_2.copy()
        return([],syndict)
    if(len(disease_STR) == 0):
        syndict = syndict_2.copy()
        return([],syndict)
    if(disease_STR == ""):
        syndict = syndict_2.copy()
        return([],syndict)
    if(disease_STR in syndict_2):
        if(syndict_2[disease_STR] == ["none"]):
            syndict = syndict_2.copy()
            return([],syndict)
        else:
            syndict = syndict_2.copy()
            return(syndict_2[disease_STR],syndict)
    #auth = get_ticket("TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas")
    auth = get_ticket("TGT-1324012-CDXW7ErncSwjrLVeYAKog5b2oMJUZtkPV3DXMYNer9l1kNs5sQ-cas")
    #auth = "ST-2227617-atRsjZqEFS6te6eP1Vme-cas"
    disease = urllib.parse.quote_plus(disease_STR)
    #disease = tmp.geturl()
    #print(disease)
    #url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=" + disease + "&ticket=" + auth
    #url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=dilated+cardiomyopathy&ticket=ST-461542-ea3bE23tuEfHThSfcfZo-cas"
    url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=" + disease + "&ticket=" + auth
    #url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=dilated+cardiomyopathy"+ "&ticket=" + auth
    #print(url)
    request = requests.get(url)
    auth = get_ticket("TGT-1324012-CDXW7ErncSwjrLVeYAKog5b2oMJUZtkPV3DXMYNer9l1kNs5sQ-cas")
    url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=" + disease + "&ticket=" + auth
    #url = "https://uts-ws.nlm.nih.gov/rest/search/current?string=dilated+cardiomyopathy"+ "&ticket=" + auth
    if not(request.status_code == 200):
        syndict_2[disease_STR] = ["none"]
        syndict = syndict_2.copy()
        return([],syndict)
    try:
        response = urllib.request.urlopen(Request(url, headers={'User-Agent': 'Mozilla/5.0'}))
        webContent = str(response.read())
        #print(webContent)
    except:
        #print("urlopen not run")
        syndict_2[disease_STR] = ["none"]
        syndict = syndict_2.copy()
        return([],syndict)
    #req = Request('http://www.cmegroup.com/trading/products/#sortField=oi&sortAsc=false&venues=3&page=1&cleared=1&group=1', headers={'User-Agent': 'Mozilla/5.0'})
    #webpage = urlopen(req).read()
    if("403" in str(response.getcode())):
        syndict_2[disease_STR] = ["none"]
        syndict = syndict_2.copy()
        return([],syndict)
    #table_file = open("VEP_test/test_syn.json")
    #webContent=table_file.read()
    if not("\"results\":[" in webContent):
        syndict_2[disease_STR] = ["none"]
        syndict = syndict_2.copy()
        return([],syndict)
    result_tmp = webContent.split("\"results\":[")[1]
    result = result_tmp.split("]")[0]
    result_list = result.split("}")
    ret = []
    if not(disease_STR in syndict_2):
        syndict_2[disease_STR] = []
    ctr = 0
    for i in result_list:
        if not("\"name\":\"" in i):
            continue
        resultSplit = i.split("\"name\":\"")[1]
        resultStr = resultSplit.split("\"")[0]
        if not(resultStr.lower() == ""):
            if not(resultStr.lower() in syndict_2):
                syndict_2[resultStr.lower()] = []
            #print("adding keys")
            #print(resultStr.lower())
            ret.append(resultStr.lower())
            syndict_2[disease_STR].append(resultStr.lower())
            syndict_2[resultStr.lower()].append(disease_STR)
            ctr = ctr + len(resultStr.lower())
            ctr = ctr + 1
    #print("returning")
    #print("disease name: ")
    #print(disease_STR)
    #print("returned: ")
    #print(ret)
    #print(str(ctr))
    #print(str(len(ret)))
    if(len(ret) == 0):
        syndict_2[disease_STR] = ["none"]
        #fh.write(disease_STR + "\tnone")
    #fh.close()
    syndict = syndict_2.copy()
    return(ret,syndict)

def get_inheritance_mode(gene_name):
    #print(gene_name)
    return("not given")
    url = "https://www.omim.org/search?index=entry&search=" + str(gene_name) + "&retrieve=geneMap&genemap_exists=true&format=tsv"
    response = urllib.request.urlopen(url)
    webContent = str(response.read())
    #print(webContent)
    #table_file = open("VEP_test/test_omim.tsv")
    #webContent=table_file.read()
    lines =webContent.split("\n")
    #print(str(len(lines)))
    if(len(lines) < 2):
        return("not given")
    has_inheritance = "false"
    contains_data = "false"
    index_of_inheritance = 0
    index_of_gene_name = 0
    for i in range(0,len(lines)):
        if(contains_data == "true" and has_inheritance == "true"):
            if(":" in lines[i]):
                lineSplit = lines[i].split("\t")
                if(gene_name in str(lineSplit[index_of_gene_name])):
                    #print(lineSplit)
                    inh = str(lineSplit[index_of_inheritance])
                    return(inh)
        if("Cytogenetic location" in lines[i]):
            #print("Cytogenetic location:")
            #print(lines[i])
            contains_data = "true"
            lineSplit = lines[i].split("\t")
            #print(lineSplit)
            for k in range(0,len(lineSplit)):
                if("Inheritance" in lineSplit[k]):
                    index_of_inheritance = k
                    has_inheritance = "true"
                #if(lineSplit[k] == "Gene/Locus"):
                if("Gene/Locus" in lineSplit[k]):
                    index_of_gene_name = k
                    #print(lineSplit)
    return("not given")

def read_pop_dict(path):
    pop_file = open(path)
    htmlstr=pop_file.read()
    lines= htmlstr.split("\n")
    ret_dict={}
    for line in lines:
        lineSplit=line.split()
        if(len(lineSplit)<2):
            continue
        if not(lineSplit[0] in ret_dict):
            ret_dict[lineSplit[0]] = lineSplit[1]
    return(ret_dict)

def read_gene_associations_2(path):
    df = pd.read_csv(path,sep=',')
    ret_dict_1 = {}
    ret_dict_2 = {}
    for col in df.columns:
        curr_col = list(df[col])
        curr_col = [i for i in curr_col if not(str(i) == "nan")]
        for i in range(0,len(curr_col)):
            if(str(curr_col[i]) in ret_dict_1):
                ret_dict_1[str(curr_col[i])].append(str(col))
            else:
                ret_dict_1[str(curr_col[i])] = []
                ret_dict_1[str(curr_col[i])].append(str(col))
            if(str(col) in ret_dict_2):
                ret_dict_2[str(col)].append(str(curr_col[i]))
            else:
                ret_dict_2[str(col)] = []
                ret_dict_2[str(col)].append(str(curr_col[i]))
    #print(ret_dict_1)
    return([ret_dict_1,ret_dict_2])

def get_hgvs_list(rsid):
    url = "https://mutalyzer.nl/snp-converter?rs_id=" + rsid
    url = "https://mutalyzer.nl/json/getdbSNPDescriptions?rs_id=" + rsid
    response = urllib.request.urlopen(url)
    webContent = response.read()
    hgvs_list = str(webContent).replace("b'","").replace("[","").replace("]","").replace('"',"").split(",")
    ret = [i.replace("'","") for i in hgvs_list]
    return ret
    #temp1 = str(webContent).split("<h4>HGVS descriptions</h4>")[1]
    #temp2 = temp1.split("</div>")[0]
    #hgvs_lines_tmp = temp2.replace("</code></p>","").replace("\\n","").replace(" ","").replace("&gt;",">").split("<p><code>")
    #hgvs_lines = [line for line in hgvs_lines_tmp if(len(line) > 1)]
    #print(hgvs_lines_tmp)
    #for line in hgvs_lines:

def get_inheritance_mode_2(hgnc_id):
    #print(hgnc_id)
    url = "https://search.clinicalgenome.org/kb/genes/HGNC:" + str(hgnc_id)
    response = urllib.request.urlopen(url)
    webContent = str(response.read())
    dfs = pd.read_html(webContent)
    if("ClinGen has not yet curated" in webContent):
        return("not given")
    elif(len(dfs) < 2):
        return("not given")
    df = dfs[2]
    #print(len(df))
    #print(list(df.columns))
    if not('MOI' in list(df.columns)):
        return("not given")
    inh_mode = df.loc[0,'MOI'].replace("\\t","").replace("\\n","")
    return(inh_mode)

def get_inheritance_mode_3(gene_name,path,disease_curr,dis_check):
    dis_syn_list = [disease_curr]
    dis_syn_list_tmp = get_synonyms(disease_curr,"TGT-1046014-OgTp1p3CqoWt5txPDftCGsUtKbGNQTdNwnOjAZ5TfKktxCxOAe-cas")
    dis_syn_list = dis_syn_list + dis_syn_list_tmp
    df = pd.read_csv(path,sep=',')
    df = df.set_index("Gene")
    if(gene_name in df.index):
        inh = df.loc[gene_name,"MOI"]
        dis = df.loc[gene_name,"Disease"]
        if(dis_check == "false"):
            if(isinstance(dis, pd.core.series.Series) and isinstance(inh, pd.core.series.Series)):
                
                inh = ";".join(inh.tolist())
            else:
                inh = str(inh)
        elif(isinstance(dis, pd.core.series.Series) and isinstance(inh, pd.core.series.Series)):
            if(len(dis.tolist()) > 1):
                if(any(dis_name in dis_syn_list for dis_name in dis.tolist())):
                    disease_found_name = next(i for i in dis_syn_list if i in dis.tolist())
                    index_curr = dis.tolist().index(disease_found_name)
                    inh = inh.tolist()[index_curr]
                else:
                    inh = "not given"
            else:
                if(dis.tolist()[0] in dis_syn_list):
                    inh=inh.tolist()[0]
                elif(dis_check == "true"):
                    inh = "not given"
                else:
                    inh=inh.tolist()[0]
        else:
            if(dis in dis_syn_list):
                inh = str(inh)
            elif(dis_check == "true"):
                inh = "not given"
            else:
                inh = str(inh)
                    
    else:
        inh = "not given"
    return(inh)

def read_gene_associations(path):
    df = pd.read_csv(path,sep=',')
    #df = df.set_index("Gene")
    ret_dict_1 = {}
    ret_dict_2 = {}
    #print(df)
    #print(list(df.columns))
    df = df.drop(df.index[[0]])
    for i in list(df.index):
        gene = str(df.loc[i,"Gene"])
        #print(gene)
        disease_curr = str(df.loc[i,"Disease"])
        #print(disease_curr)
        if not(gene in ret_dict_1):
            ret_dict_1[gene] = []
        ret_dict_1[gene].append(disease_curr)
        if not(disease_curr in ret_dict_2):
            ret_dict_2[disease_curr] = []
        ret_dict_2[disease_curr].append(gene)
    return([ret_dict_1,ret_dict_2])
            
    
def get_link(hgnc_id,disease):
    url = "https://search.clinicalgenome.org/kb/genes/HGNC:" + str(hgnc_id)
    response = urllib.request.urlopen(url)
    webContent = str(response.read())
    if("ClinGen has not yet curated" in webContent):
        return("N/A")
    if not("Gene-Disease Validity</h3>" in webContent):
        return("N/A")
    webContentSplit = webContent.split("Gene-Disease Validity</h3>")[1]
    dfs = pd.read_html(webContentSplit)
    if(len(dfs) < 1):
        return("N/A")
    df = dfs[0]
    nbr_col = 0
    #print(list(df.columns))
    #print(df.loc[0,'Report & Date'])
    if("Report & Date" in list(df.columns)):
        nbr_col = list(df.columns).index("Report & Date")
        table_html = webContentSplit.split("</table>")[0]
        table_html_2 = table_html.split("<tbody class=")[1]
        table_rows = table_html_2.split("<td class=")
        link = table_rows[4].split("href=\"")[1].split("\">")[0]
    url = "https://search.clinicalgenome.org/" + str(link)
    return(link)
        #return(table_row_with_link)
        #link = df.loc[0,'Report &amp; Date']
        #return(link)
    
    
def get_publications(rsid):
    ret = {}
    ids_STR = ""
    titles_and_ids_STR = ""
    url = "https://www.ncbi.nlm.nih.gov/snp/" + rsid
    response = urllib.request.urlopen(url)
    request = requests.get(url)
    if not(request.status_code == 200):
        return("not given")
    webContent = str(response.read())
    publ_div_tmp = webContent.split("<div id=\"publication\">")
    if(len(publ_div_tmp) < 2):
        return([{"-9999":"-9999"},"-9999","-9999"])
    publ_div = publ_div_tmp[1].split("</div>")[0]
    #print(publ_div)
    dfs = pd.read_html(publ_div)
    df = dfs[0]
    #print(df)
    #print(list(df.columns))
    ids = []
    for i in range(0,len(df.index)):
        pmid = df.loc[i,'PMID']
        title = df.loc[i,'Title']
        pmid = pmid.rstrip()
        pmid = pmid.split("\\")[0]
        #print(pmid)
        #if("href=\"/pubmed/" in pmid):
        #    pubmed_tmp = pmid.split("href=\"/pubmed/")[1].split("\"")[0]
        #    pmid_str = str(pubmed_tmp)
        ret[pmid] = str(title)
        ids.append(pmid)
    ids_STR = str( ";".join(ids))
    titles_and_ids_STR = ";".join([str(ret[key]) + "(PMID: " + str(key) + ")" for key in ret])
    return([ret,ids_STR,titles_and_ids_STR])
       


    
    
    
    
def get_correct_hgvs_link(hgvs_list,gene_name,prot_change):
    for i in hgvs_list:
        #print(i)
        temp1 = i.split(":")[0]
        temp2 = i.split(":")[1]
        querystr= temp1 + "(" + gene_name + ")%3A" + temp2 + "%20(" + prot_change + ")"
        querystr = querystr.replace(" ","")
        #print(querystr)
        query = "https://clinvarminer.genetics.utah.edu/submissions-by-variant/" + querystr
        #print(query)
        #print(query)
        request = requests.get(query)
        #request = requests.head(query)
        if request.status_code == 200:
            return(query)
    return(["",""])
        
def get_correct_hgvs(hgvs_list,gene_name,prot_change):
    for i in hgvs_list:
        #print(i)
        temp1 = i.split(":")[0]
        temp2 = i.split(":")[1]
        querystr= temp1 + "(" + gene_name + ")%3A" + temp2 + "%20(" + prot_change + ")"
        querystr = querystr.replace(" ","")
        #print(querystr)
        query = "https://clinvarminer.genetics.utah.edu/submissions-by-variant/" + querystr
        #print(query)
        #print(query)
        request = requests.get(query)
        #request = requests.head(query)
        if request.status_code == 200:
            return([i,query])
    return(["",""])
            #os.system("python3 extract_table.py \"" + query + "\"")
            #break

#link = get_correct_hgvs_link(get_hgvs_list(snpid),gene_name,prot_change)
#print(link)
#ret = extract(link)
#for key in ret:
#    print(str(key) + ":\t" + str(ret[key]) + "\n")

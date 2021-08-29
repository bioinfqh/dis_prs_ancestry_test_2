
import sys
import numpy as np
import os
import pandas as pd
from scipy import stats
from sklearn.metrics import confusion_matrix

#from Utils.postprocess import get_samples_from_msp_df
from Utils.visualization import plot_cm, plot_chm


xgmix_path = "/XGMix_master/XGMIX.py"

def run_calc(path_to_infiles,path_to_models,output_basename):
    outfiles = []
    for i in range(21,23):
        curr_chr = str(i)
        curr_model = path_to_models + "/chm_" + curr_chr + ".pkl.gz"
        curr_vcf = path_to_infiles + "/chr_" + curr_chr + ".vcf.gz"
        #print("python3 XGMIX.py " + curr_vcf + " demo_data/allchrs.b37.gmap " +  output_basename + "_chr_" + curr_chr + " " + curr_chr + " False " + curr_model)
        if(i == 21):
            os.system("python3 " + xgmix_path + " " + curr_vcf + " /XGMix-master/demo_data/allchrs.b37.gmap " +  output_basename + "_chr_" + curr_chr + " " + curr_chr + " False " + curr_model)
        #print(output_basename + "_chr_" + curr_chr)
        if(os.path.isdir(output_basename + "_chr_" + curr_chr)):
            #print(" " + output_basename + "_chr_" + curr_chr + "/demo_chr_" + curr_chr + ".msp.tsv " + output_basename)
            #print("cp " + output_basename + "_chr_" + curr_chr + "/demo_chr_" + curr_chr + ".msp.tsv demo_data/demo_chr_" + curr_chr + ".msp.tsv")
            os.system("cp " + output_basename + "_chr_" + curr_chr + "/demo_chr_" + curr_chr + ".msp.tsv /XGMix-master/demo_data/demo_chr_" + curr_chr + ".msp.tsv")
        curr_output = output_basename + "_chr_" + curr_chr
        outfiles.append(curr_output)
    return(outfiles)

def read_results_and_plot(sample_id,output_file):
    msp_df = pd.read_csv(output_file, sep="\t", skiprows=[0])
    #print(msp_df)
    #msp_df.loc[:,'#chm'] = 1
    #print(msp_df)
    #img_name = "demo_data/imgs/chm_img_" + sample_id
    img_name = "/XGMix-master/demo_data/imgs/chm_img"
    plot_chm(sample_id, msp_df, img_name)
    os.system('mv /XGMix-master/demo_data/chm_img.svg /XGMix-master/demo_data/chm_img_' + sample_id + '.svg')
    os.system('mv /XGMix-master/demo_data/chm_img.png /XGMix-master/demo_data/chm_img_' + sample_id + '.png')
    return(img_name)


def plot_all(sample_id,output_basename):
    #msp_df = pd.read_csv((output_basename + "_chr_1.msp.tsv"), sep="\t", skiprows=[0])
    msp_df = pd.read_csv((output_basename + "_chr_21.msp.tsv"), sep="\t", skiprows=[0])
    #print(output_basename + "_chr_1.msp.tsv")
    #msp_df = pd.read_csv((output_basename + ".msp.tsv"), sep="\t", skiprows=[0])
    for i in range(22,23):
        curr_chr = str(i)
        #print(output_basename + "_chr_" + curr_chr + ".msp.tsv")
        curr_msp_df = pd.read_csv((output_basename + "_chr_" + curr_chr + ".msp.tsv"), sep="\t", skiprows=[0])
        #print(curr_msp_df)
        msp_df_new = pd.concat([msp_df, curr_msp_df], ignore_index=True)
        msp_df = msp_df_new.copy()
    #img_name = "demo_data/imgs/chm_img_" + sample_id
    img_name = "/XGMix_master/demo_data/imgs/chm_img"
    #print(msp_df)
    plot_chm(sample_id, msp_df, img_name)
    os.system('mv chm_img.svg chm_img_' + sample_id + '.svg')
    os.system('mv chm_img.png chm_img_' + sample_id + '.png')
    return(["chm_img_" + sample_id + ".svg","chm_img_" + sample_id + ".png"])

def make_html(img_paths):
    html_file = open("local_ancestry_to_fill.html")
    htmlstr=html_file.read()
    html_file.close()
    if(len(img_paths) < 1):
        return(htmlstr)
    htmlstr_temp = htmlstr.replace("chrom_anc_img",img_paths[1])
    htmlstr = htmlstr_temp
    return(htmlstr)


sample_id = "NA21141"
#read_results_and_plot(sample_id,"demo_data/demo.msp.tsv")

def run_and_plot(path_to_infiles,sample_id,path_to_models,output_basename):
    outfiles = run_calc(path_to_infiles,path_to_models,output_basename)
    plot_all(sample_id,output_basename)

def run_and_html(path_to_infiles,sample_id,path_to_models,output_basename,customer_id):
    outfiles = run_calc(path_to_infiles,path_to_models,output_basename)
    img_paths = plot_all(sample_id,output_basename)
    #html_str = make_html(img_paths)
    #outfile = "report_chrom_anc_" + customer_id + ".html"
    os.system('mv /XGMix-master/demo_data/chm_img_' + sample_id + '.png /scripts/chm_img_' + customer_id + '.png')
    #fh=open(outfile,'w')
    #fh.write(html_str)
    #fh.close()
    return("")

if(len(sys.argv) > 5):
    vcf_dir=str(sys.argv[1])
    sample_id=str(sys.argv[2])
    model_file=str(sys.argv[3])
    output_prefix_path=str(sys.argv[4])
    customer_id=str(sys.argv[5])
    run_and_html(vcf_dir,sample_id,model_file,output_prefix_path,customer_id)
    print("/scripts/chm_img_" + customer_id + ".png")
else:
    run_and_html("/XGMIX_infiles",sample_id,"/XGMIX_model_files","/XGMix-master/demo_data/demo","testuser")
#run_and_plot("/media/quirin/INTENSO/XGMIX_infiles",sample_id,"/media/quirin/INTENSO/XGMIX_model_files","demo_data/demo")

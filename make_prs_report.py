
#from plot_curve import make_risk_plot
#from plot_curve import make_risk_plot_2
import sys
#import pdfkit 
#import imgkit
#import os
import json
#from risk_calc_2 import get_absolute_risk_and_percentile

colorscale_width = 1100
sep = "<hr>"
#disease_name_input = sys.argv[1]
risk_limit = 0.9
params_from_file = "true"
part1 = "part1.html"
part2 = "part2.html"

def export_prs_as_json(patient_id,disease_names,risks,percs):
    dict_for_ret = {}
    for i in range(0,len(disease_names)):
        new_dict = {}
        new_dict["risk_score"] = str(risks[i])
        new_dict["percentile"] = str(percs[i])
        dict_for_ret[disease_names[i]] = new_dict
    json_str = json.dumps(dict_for_ret)
    return(json_str)

        

def make_html_easy(html_path,outfile,results_path,patient_id,risk_threshold,is_abs_risk_file,customer_id):
    disease_names = []
    percs = []
    risks = []
    distrs = []
    ## read result file
    for line in resultlist:
        lineSplit = line.split()
        if(len(lineSplit) < 4):
            continue
        else:
            disease_names.append(lineSplit[0])
            risks.append(float(lineSplit[1]))
            percs.append(float(lineSplit[2]))
            distrs.append(str(lineSplit[3]))
    json_str = export_prs_as_json(patient_id,disease_names,risks,percs)
    #print(json_str)
    endstrs = []
    #print(header_str)
    fh=open("prs_" + patient_id + ".json",'w')
    fh.write(json_str)
    fh.close()
    return("prs_" + patient_id + ".json")        
    
    
if(params_from_file == "true"):
    #imgkit.from_file("gauge.html","new_scale.jpg")
    #make_scale_image("gauge.html","new_scale.jpg")
    #pdfkit.from_file("gauge.html","new_scale.pdf")
    print(sys.argv[1])
    ## this is the main mode that is actually used
    if(sys.argv[1] == "easy"):
        result_file=str(sys.argv[2])
        patient_id=str(sys.argv[3])
        output_file = str(sys.argv[4])
        is_abs_risk = str(sys.argv[5])
        #outfile_curr = outfile + "_prs.html"
        print(output_file)
        patient_id_new = output_file.replace("prs_report_","").replace(".pdf","")
        make_html_easy("PRS_to_fill.html","html_page.html",result_file,patient_id,float("0.7"),is_abs_risk,patient_id_new)
        #make_pdf_from_str(make_html_easy("PRS_to_fill.html","html_page.html",result_file,patient_id,float("0.7"),is_abs_risk,patient_id_new),output_file)
        #make_pdf_from_str(make_html_easy("PRS_to_fill.html","none",result_file,patient_id,float("0.7"),is_abs_risk),output_file)

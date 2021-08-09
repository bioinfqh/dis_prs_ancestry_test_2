import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import sys


simulate="false"
col_nbr_with_score=2
patient_id_temp="HG00096"
col_with_ids=1
#patient_risk = 0.0
plot_risk="false"

def make_risk_plot(path,patient_id,outfile):
    patient_risk = 0.0
    fig, ax = plt.subplots(figsize=(12, 6))
    #x_values = np.linspace(-3, 3, 120)
    if(simulate == "true"):
        x_axis = np.arange(0, 100, 0.001)
        ax.plot(x_axis, norm.pdf(x_axis,0,2))
        patient_risk_percentile=0.5
    else:
        results = open(path)
        htmlstr=results.read()
        lines= htmlstr.split("\n")
        y = []
        for i in range(1,len(lines)):
            if(len(lines[i].split()) < col_with_ids):
                continue
            #print(lines[i].split())
            if(lines[i].split()[col_with_ids-1] == patient_id and not(len(lines[i].split()) < col_nbr_with_score)):
                patient_risk = float(str(lines[i].split()[col_nbr_with_score-1]))
                print(patient_risk)
            if not(len(lines[i].split()) < col_nbr_with_score):
                y.append(float(str(lines[i].split()[col_nbr_with_score-1])))
        #print(y)
        #print(np.linspace(min(y),max(y),((max(y) - min(y)) / 100.0)))
        #print(str(min(y)))
        #print(str(max(y)))
        hist,bins=np.histogram(y,bins=np.linspace(min(y),max(y),100))
        #bins_new = [((float(bins[i]) + float(bins[i+1]))/2.0) for i in range(0,len(bins)-1)]
        bins_new = []
        y_len = len(y)
        hist_new = [(float(i) / float(y_len)) for i in hist]
        for i in range(0,len(bins)-1):
            bins_new.append((bins[i] + bins[i+1]) / 2.0)
        ax.plot(bins_new,hist_new)
        lower_risk=[i for i in y if (i < patient_risk)]
        patient_risk_percentile = float(len(lower_risk))/float(y_len)
        patient_risk_score = ((patient_risk - min(y)) / (max(y)-min(y)))
    ax.vlines(patient_risk, 0, max(hist_new), linestyles='dashed', colors='red')
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    fig.savefig(outfile)
    return(patient_risk,patient_risk_percentile,patient_risk_score)
    #patient_risk_percentile = ((patient_risk - min(y)) / (max(y)-min(y)))
    

def make_risk_plot_2(risk_dict,patient_id,outfile):
    patient_risk = 0.0
    fig, ax = plt.subplots(figsize=(12, 6))
    #x_values = np.linspace(-3, 3, 120)
    if(simulate == "true"):
        x_axis = np.arange(0, 100, 0.001)
        ax.plot(x_axis, norm.pdf(x_axis,0,2))
        patient_risk_percentile=0.5
    else:
        #results = open(path)
        #htmlstr=results.read()
        #lines= htmlstr.split("\n")
        y = []
        patient_risk = float(risk_dict[patient_id])
        y = [risk_dict[i] for i in risk_dict]
        #print(y)
        #print(np.linspace(min(y),max(y),((max(y) - min(y)) / 100.0)))
        hist,bins=np.histogram(y,bins=np.linspace(min(y),max(y),100))
        #bins_new = [((float(bins[i]) + float(bins[i+1]))/2.0) for i in range(0,len(bins)-1)]
        bins_new = []
        y_len = len(y)
        hist_new = [(float(i) / float(y_len)) for i in hist]
        for i in range(0,len(bins)-1):
            bins_new.append((bins[i] + bins[i+1]) / 2.0)
        ax.plot(bins_new,hist_new)
        lower_risk=[i for i in y if (i < patient_risk)]
        patient_risk_percentile = float(len(lower_risk))/float(y_len)
        patient_risk_score = ((patient_risk - min(y)) / (max(y)-min(y)))
    ax.vlines(patient_risk, 0, max(hist_new), linestyles='dashed', colors='red')
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    fig.savefig(outfile)
    return(patient_risk,patient_risk_percentile,patient_risk_score)
    #patient_risk_percentile = ((patient_risk - min(y)) / (max(y)-min(y)))
    
#perc = make_risk_plot(sys.argv[1],sys.argv[2],sys.argv[3])
#print(perc)

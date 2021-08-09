/*
 * relfinder_vcf.cpp
 *
 *  Created on: Nov 8, 2017
 *      Author: joepickrell
 */

#include "DataVCF.h"

void printv(){
	cout << "\nrelfinder_vcf v. 0.00\n";
	cout << "Joe Pickrell (joepickrell@gmail.com)\n\n";
}

void printopts(){
	cout << "Options:\n";
	cout << "-vcf1 [file name] a VCF file with personal allele frequencies and allele counts (from ancestry).\n";
	cout << "-vcf2 [file name] a second VCF file with personal allele frequencies and allele counts (from ancestry).\n";
	cout << "-o [string] stem for output files (relfinder)\n";
	cout << "\n";
}


int main(int argc, char *argv[]){
	printv();
    CCmdLine cmdline;
    VCF_params p;
    if (cmdline.SplitLine(argc, argv) < 1){
    		printopts();
    		exit(1);
    }
    string infile1, infile2;
    string outstem = "relfinder";
    if (cmdline.HasSwitch("-vcf1")) p.infile1 = cmdline.GetArgument("-vcf1", 0);
    else{
    		printopts();
    		exit(1);
    }
    if (cmdline.HasSwitch("-vcf2")) p.infile2= cmdline.GetArgument("-vcf2", 0);
    else{
    		printopts();
    		exit(1);
    }
    if (cmdline.HasSwitch("-o")) outstem = cmdline.GetArgument("-o", 0);

    DataVCF d(&p);

    // run the inference
    vector<pair<double, double> > rel = d.relfinder();
    string oflk = outstem+".llk";
    ofstream outf(oflk.c_str());
    for (vector<pair<double, double> >::iterator it = rel.begin(); it != rel.end(); it++){
    	outf << it->first << " "<< it->second << "\n";
    }

    // get maximum
    int maxi = 0;
    float maxlk = rel[0].second;
    for (int i = 0; i < rel.size(); i++){
    		if (rel[i].second > maxlk){
    			maxi = i;
    			maxlk = rel[i].second;
    		}
    }
    //float maxlk = rel[maxi].second;
    // get CIs
    pair<float, float> ci;
    int loi = maxi;
    int hii = maxi;
    bool lofound = false;
    bool hifound = false;
    while (loi >0 and !lofound){
    		float tmplk = rel[loi].second;
    		float lr = 2.0*(maxlk-tmplk);
    		if (lr > 3.84){
    			loi++;
    			lofound = true;
    		}
    		loi --;
    }
    while (hii< rel.size() and !hifound){
    		float tmplk = rel[hii].second;
    		float lr = 2.0*(maxlk-tmplk);
    		if (lr > 3.84){
    			hii--;
    			hifound = true;
    		}
    		hii ++;
    }
    double lr = 2.0*(maxlk - rel[0].second);
   // cout << "lr "<< lr << "\n"; cout.flush();
    double pval = 1-gsl_cdf_chisq_P (lr, 1.0);
    string ofsummary = outstem+".summary";
    ofstream outf2(ofsummary.c_str());
    outf2 <<"{\n";
    outf2 << "\t\"est\": "<< rel[maxi].first << ",\n";
    outf2 << "\t\"lo\": "<< rel[loi].first << ",\n";
    outf2 << "\t\"hi\": "<< rel[hii].first << ",\n";
    outf2 << "\t\"maxlk\": "<< rel[maxi].second <<",\n";
    outf2 << "\t\"lk0\": "<< rel[0].second << ",\n";
    outf2 << "\t\"nsnp\": "<< d.overlapindex.size() << ",\n";
    outf2 << "\t\"P\": "<< pval<< "\n";
    outf2 <<"}\n";

    return 0;

}



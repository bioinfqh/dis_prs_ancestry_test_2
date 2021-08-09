/*
 * DataVCF.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: jkpickrell
 *
 *
 *      The squarem1() and squarem2() functions are modified from https://github.com/daichaoxing/cxxSQUAREM
 *
 *
 */
#include "DataVCF.h"

DataVCF::DataVCF(){

}


DataVCF::DataVCF(VCF_params *p){
	params = p;
	data.clear();
	info.clear();
	chr2pos2index.clear();
	f1.clear();
	overlapindex.clear();

	read_vcf1(p->infile1);
	read_vcf2(p->infile2);


}


vector<pair<double, double> > DataVCF::relfinder(){
	vector<pair<double, double> > toreturn;
	vector<double> weighted_alfreq2;
	cout << "### relfinder ###" <<"\n";
	cout << "### overlap of " << overlapindex.size()<< " SNPs\n";
	cout <<"###\n";
	for (double i = 0; i <=1.001; i+=0.01){
		double llk = relfinder_llk(i);
		toreturn.push_back(make_pair(i, llk));
	}
	return toreturn;
}

double DataVCF::relfinder_llk(float lambda){
	double toreturn = -10000;
	for (int i = 0; i < overlapindex.size() ;i++){
		int index = overlapindex[i];
		pair<bool, bool> dd = data.at(index);
		assert(f1.at(index) >=0);
		bool g1 = dd.first;
		bool g2 = dd.second;
		bool eq = false;
		double f = f1.at(index);
		if (!g2) f = 1-f;
		if (g1 == g2) eq = true;
		double pmatch = lambda *(0.5 + 0.5* f) + (1-lambda) *f;
		if (pmatch < 1e-10) pmatch = 1e-10;
		if (pmatch > 1-1e-10) pmatch = 1-1e-10;
		if (eq) toreturn += log(pmatch);
		else toreturn+= log(1.0-pmatch);

	}
	return toreturn;

}




void DataVCF::read_vcf1(string infile){
	// read a VCF/BCF file
	//
	// not checking if genome build is the same as in the reference frequencies
	// checking to make sure the allele matches one of the reference alleles
	cout << "Reading "<< infile << "\n";
    int nwarning = 0;
	//using htslib
	bcf_srs_t *sr =  bcf_sr_init();
	bcf_sr_add_reader (sr, infile.c_str() );
	bcf_hdr_t *header = sr->readers->header;
	bcf1_t *line; //VCF line gets put in here
	int Nsamp=bcf_hdr_nsamples(header);
	if (Nsamp > 1){
		cerr << "ERROR: should only be one sample per file, this file has "<< Nsamp << " samples\n";
		exit(1);
	}
	int index = 0;
	while(bcf_sr_next_line (sr)) { //loop through file
	   line =  bcf_sr_get_line(sr, 0);  //read a line
	   bcf_unpack(line, BCF_UN_ALL);

	   string chr =  bcf_hdr_id2name(header, line->rid);
	   int pos = line->pos+1;

	   // ref counts
	   int nrc_arr = 0;
	   int nrc     = 0;
	   int *rc     = NULL;
	   // alt counts
	   int nac_arr = 0;
	   int nac     = 0;
	   int *ac     = NULL;

	   // altf
	   int niaf_arr = 0;
	   float niaf   = 0.0;
	   float *iaf     = NULL;
	   nrc = bcf_get_format_int32(header, line, "RC", &rc, &nrc_arr);
	   nac = bcf_get_format_int32(header, line, "AC", &ac, &nac_arr);
	   niaf = bcf_get_format_float(header, line, "IAF", &iaf, &niaf_arr);
      // cout << "ok\n"; cout.flush();
     //  cout << rc[0] << " "<< ac[0] << " "<< iaf[0] << " "<< chr << " "<< pos <<"\n";
       info.push_back(make_pair(chr, pos));
       if (chr2pos2index.find(chr) == chr2pos2index.end()){
    	   	   map<int, int> tmp;
    	   	   tmp[pos]= index;
    	   	   chr2pos2index.insert(make_pair(chr, tmp));
       }
       else{
    	   	   chr2pos2index[chr].insert(make_pair(pos, index));
       }
       bool dac = false;
       if (ac[0] >0) dac = true;
       bool dac2;
       data.push_back(make_pair(dac, dac2));
       f1.push_back(iaf[0]);
       index++;
	}
	free(sr);
	free(line);
	free(header);
}

void DataVCF::read_vcf2(string infile){
	// read a VCF/BCF file
	//
	// not checking if genome build is the same as in the reference frequencies
	// checking to make sure the allele matches one of the reference alleles
	cout << "Reading "<< infile << "\n";
    int nwarning = 0;
	//using htslib
	bcf_srs_t *sr =  bcf_sr_init();
	bcf_sr_add_reader (sr, infile.c_str() );
	bcf_hdr_t *header = sr->readers->header;
	bcf1_t *line; //VCF line gets put in here
	int Nsamp=bcf_hdr_nsamples(header);
	if (Nsamp > 1){
		cerr << "ERROR: should only be one sample per file, this file has "<< Nsamp << " samples\n";
		exit(1);
	}
	while(bcf_sr_next_line (sr)) { //loop through file
	   line =  bcf_sr_get_line(sr, 0);  //read a line
	   bcf_unpack(line, BCF_UN_ALL);

	   string chr =  bcf_hdr_id2name(header, line->rid);
	   int pos = line->pos+1;

	   // ref counts
	   int nrc_arr = 0;
	   int nrc     = 0;
	   int *rc     = NULL;
	   // alt counts
	   int nac_arr = 0;
	   int nac     = 0;
	   int *ac     = NULL;

	   // altf
	   int niaf_arr = 0;
	   float niaf   = 0.0;
	   float *iaf     = NULL;
	   nrc = bcf_get_format_int32(header, line, "RC", &rc, &nrc_arr);
	   nac = bcf_get_format_int32(header, line, "AC", &ac, &nac_arr);
	   niaf = bcf_get_format_float(header, line, "IAF", &iaf, &niaf_arr);

       //cout << rc[0] << " "<< ac[0] << " "<< iaf[0] << " "<< chr << " "<< pos <<"\n";
       if (chr2pos2index.find(chr) != chr2pos2index.end() and chr2pos2index[chr].find(pos) != chr2pos2index[chr].end()){
    	   	   int index = chr2pos2index[chr][pos];
    	       bool dac = false;
    	       if (ac[0] >0) dac = true;
    	       data.at(index).second = dac;
    	       overlapindex.push_back(index);
       }

	}
	free(sr);
	free(line);
	free(header);
}


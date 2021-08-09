/*
 * DataVCF.h
 *
 *  Created on: Sep 20, 2014
 *      Author: jkpickrell
 */

#ifndef DATAVCF_H_
#define DATAVCF_H_
#include "VCF_params.h"

class DataVCF{
public:
	DataVCF();
	DataVCF(VCF_params*);
	VCF_params *params;

	vector<pair<string, int> > info; //vector of [chr, pos]
	map<string, map<int, int> > chr2pos2index;
	vector<pair<bool, bool> > data;
	vector<float> f1;
	vector<int> overlapindex;
	void read_vcf1(string);
	void read_vcf2(string);

	// relfinder functions
	//vector<bool> relfinder_overlap(DataVCF *);
	vector<pair<double, double> > relfinder();
	double relfinder_llk(float);

};



#endif /* DATA23_H_ */

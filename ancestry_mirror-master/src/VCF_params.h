/*
 * Ancestry_params.h
 *
 *  Created on: Sep 16, 2014
 *      Author: jkpickrell
 */

#ifndef ANCESTRY_PARAMS_H_
#define ANCESTRY_PARAMS_H_
#include <string>
#include <map>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <stack>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include "gzstream.h"
#include "CmdLine.h"
#include <sys/stat.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <algorithm>
#include <htslib/synced_bcf_reader.h>
#include <htslib/hts.h>

using std::string;
using std::vector;
using std::list;
using std::stack;
using std::map;
using std::set;
using std::multiset;
using std::cout;
using std::cin;
using std::endl;
using std::ostream;
using std::ofstream;
using std::stringstream;
using std::pair;
using std::iterator;
using std::pair;
using std::make_pair;
using std::fstream;
using std::ifstream;
using boost::math::binomial_coefficient;



class VCF_params{
public:
	VCF_params();
	string infile1, infile2;
	string outstem;
	void print();

};


#endif /* ANCESTRY_PARAMS_H_ */

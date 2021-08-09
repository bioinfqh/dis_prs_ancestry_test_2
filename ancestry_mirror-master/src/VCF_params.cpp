/*
 *VCF_params.cpp
 *
 *  Created on: nov 8 2017
 *      Author: jkpickrell
 */

#include "VCF_params.h"

VCF_params::VCF_params(){
	outstem = "relfinder";
}


void VCF_params::print(){

	cout << ":: output stem: "<< outstem << "\n";

}

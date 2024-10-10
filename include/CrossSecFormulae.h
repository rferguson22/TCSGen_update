#ifndef CROSSSECFORMULAE_H
#define CROSSSECFORMULAE_H

#include "CrossSecMaster.h"
#include <memory>
#include <vector>
#include <string>

class BH : public cross_section {
public:
	BH(shared_ptr<cross_section> target){}

	double c_sec(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th) override{
		double f1p = instance->f1(a_t);
		double f2p = instance->f2(a_t);
		double m = instance->mass;		

		return ttcscrs.Eval_BH(a_s,a_Q2,a_t,a_weight,a_phi,a_th,f1p,f2p,m);
	}
};

vector<string> get_valid_cross_sec(){
	return {"BH"};
}

#endif //CROSSSECFORMULAE_H

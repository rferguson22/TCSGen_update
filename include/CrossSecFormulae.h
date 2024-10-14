#ifndef CROSSSECFORMULAE_H
#define CROSSSECFORMULAE_H

#include "CrossSecMaster.h"
#include <memory>
#include <vector>
#include <string>

class BH : public cross_section {
private:
	std::shared_ptr<cross_section> target;
public:
	BH(std::shared_ptr<cross_section> target):target(target){}

	double c_sec(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th,double a_sc_D){
		double f1p = target->f1(a_t);
		double f2p = target->f2(a_t);
		double m = target->mass;		

		return ttcscrs.Eval_BH(a_s,a_Q2,a_t,a_weight,a_phi,a_th,f1p,f2p,m);
	}

	using cross_section::c_sec;
};

class INT : public cross_section {
private:
        std::shared_ptr<cross_section> target;
public:
        INT(std::shared_ptr<cross_section> target):target(target){}

        double c_sec(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th,double a_sc_D){
                double f1p = target->f1(a_t);
                double f2p = target->f2(a_t);
                double m = target->mass;

                return ttcscrs.Eval_INT(a_s,a_Q2,a_t,a_weight,a_phi,a_th,a_sc_D,f1p,f2p,m);
        }

        using cross_section::c_sec;
};


std::vector<std::string> get_valid_cross_sec(){
	return {"BH","INT"};
};

#endif //CROSSSECFORMULAE_H

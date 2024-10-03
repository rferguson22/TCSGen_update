#ifndef TARGETS_H
#define TARGETS_H

#include <TTCSCrs.h>

class cross_sections {
public:
	int PID;
	double mass;
	double amm;
	double c_sec(){
		return 0;
		};
	double f1(){
		return 0;
		};
	double f2(){
		return 0;
		};
};

class p_BH:public cross_sections {
public:
	TTCSCrs ttcscrs;

	p_BH();

	double c_sec(double a_s,double a_Q2, double a_t, double a_weight, double a_phi, double a_th);

	double f1(double t);
};


#endif 

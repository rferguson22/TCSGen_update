#ifndef TARGETS_H
#define TARGETS_H

#include <TTCSCrs.h>

class cross_sections {
public:
	TTCSCrs ttcscrs;

	int PID;
	double mass;
	double mag_mom;
	virtual double c_sec(){return 0;};
	virtual double c_sec(double a_s,double a_Q2, double a_t, double a_weight, double a_phi, double a_th){return 0;};

	virtual double f1(){return 0;};
	virtual double f2(){return 0;};
};

class p_BH:public cross_sections {
public:
	p_BH();

	double c_sec(double a_s,double a_Q2, double a_t, double a_weight, double a_phi, double a_th) override;

	double f1(double t);

	double f2(double t);
};

class n_BH:public cross_sections {
public:
	n_BH();

	double c_sec(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th) override;


	double f1(double t);

	double f2(double t);
};








#endif // TARGETS_H 

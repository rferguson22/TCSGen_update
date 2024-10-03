#include <TF1.h>
#include <TF2.h>
#include <Targets.h>
#include "TTCSCrs.h"

//Creation of proton target class
class p_tar : public cross_sections {
public:
	p_tar();
	double f1(double t);
};

p_tar::p_tar(){
	mass = 0.9383;
	PID = 2212;
}

double p_tar::f1(double t){
	const double m_p=mass;
	return (1/((1-t/0.71)*(1-t/0.71)))*
		(1/(1-t/(4*m_p*m_p)))*
		(1-2.79*t/(4*m_p*m_p));
}


//Implementation of bette-heidler cross section for proton target
p_BH::p_BH(){
	p_tar tar;
	mass=tar.mass;
	PID=tar.PID;
}

double p_BH::c_sec(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th){
	return ttcscrs.Eval_BH(a_s,a_Q2,a_t,a_weight,a_phi,a_th);
}


























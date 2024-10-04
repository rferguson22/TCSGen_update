#include <TF1.h>
#include <TF2.h>
#include "Targets.h"
#include "TTCSCrs.h"

//Creation of proton target class
class p_tar : public cross_sections {
public:
	p_tar();
	double f1(double t);
	double f2(double t);
};

p_tar::p_tar(){
	mass = 0.9383;
	PID = 2212;
	mag_mom = 2.793;
}

double p_tar::f1(double t){
	return (1/((1-t/0.71)*(1-t/0.71)))*
		(1/(1-t/(4*mass*mass)))*
		(1-mag_mom*t/(4*mass*mass));
}

double p_tar::f2(double t){
	return (1/((1-t/0.71)*(1-t/0.71)))*
		(1/(1-t/(4*mass*mass)))*
		(mag_mom-1);
}


	

//Implementation of bette-heidler cross section for proton target
p_BH::p_BH(){
	p_tar tar;
	mass=tar.mass;
	PID=tar.PID;
}

double p_BH::f1(double t){return p_tar().f1(t);}
double p_BH::f2(double t){return p_tar().f2(t);}

double p_BH::c_sec(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th){
	double f1p=this->f1(a_t);
	double f2p=this->f2(a_t);

	return ttcscrs.Eval_BH(a_s,a_Q2,a_t,a_weight,a_phi,a_th,f1p,f2p,mass);
}


























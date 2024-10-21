#ifndef CROSSSECMASTER_H
#define CROSSSECMASTER_H

#include <TTCSCrs.h>
#include <string>
#include <memory>
#include <vector>
#include <iostream>

class cross_section {
public:
        int PID;
        double mass;
        double mag_mom;

	TTCSCrs ttcscrs;
	
	virtual ~cross_section()=default;

        virtual double c_sec(){return 0;};
	virtual double c_sec(double a_s, double a_Q2, double a_t,double a_weight, double a_phi, double a_th,double a_sc_D ){return 0;};

        virtual double f1(double t){return 0;};
        virtual double f2(double t){return 0;};
	virtual std::string get_type() const=0;

	virtual int get_PID() const {return PID;};
        virtual double get_mass() const {return mass;};
        virtual double get_mag_mom() const {return mag_mom;};
};
	

#endif //CROSSSECMASTER_H

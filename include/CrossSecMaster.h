#ifndef CROSSSECMASTER_H
#define CROSSSECMASTER_H

#include <TTCSCrs.h>
#include <string>
#include <memory>
#include <vector>
#include <iostream>

class cross_section {
protected:
        int PID;
        double mass;
        double mag_mom;
public:
	TTCSCrs ttcscrs;
	
	virtual ~cross_section()=default;

        virtual double c_sec(){return 0;};
        virtual double f1(){return 0;};
        virtual double f2(){return 0;};
	virtual std::string get_type() const =0;

};


#endif //CROSSSECMASTER_H

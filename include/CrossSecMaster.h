#ifndef CROSSSECMASTER_H
#define CROSSSECMASTER_H

#include <TTCSCrs.h>
#include <string>

class cross_sections {
public:
        TTCSCrs ttcscrs;

        int PID;
        double mass;
        double mag_mom;
        virtual double c_sec(){return 0;};
        virtual double f1(){return 0;};
        virtual double f2(){return 0;};
	vitrual std::string className() const=0;
};


#endif //CROSSSECMASTER_H

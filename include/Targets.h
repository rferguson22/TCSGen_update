#ifndef TARGETS_H
#define TARGETS_H

#include "CrossSecMaster.h"  

class p_tar : public cross_section {
public:
	p_tar(){
		mass = 0.9383;
		PID = 2212;
		mag_mom = 2.793;
	}
	
protected:
	double mass;
	int PID;
	double mag_mom;
};

class n_tar : public cross_section {
public:
	n_tar(){
                mass = 0.939565;
                PID = 2112;
                mag_mom = -1.913;
        }

protected:
        double mass;
        int PID;
        double mag_mom;
};


vector<string> valid_targets = {"p_tar","n_tar"};

#endif // TARGETS_H

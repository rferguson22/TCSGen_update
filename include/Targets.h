#ifndef TARGETS_H
#define TARGETS_H

#include "CrossSecMaster.h"  
#include <vector>
#include <string>


class p_tar : public cross_section {
public:
	p_tar(){
		mass =0.938272;
		PID = 2212;
		mag_mom = 2.793;
	}
};

class n_tar : public cross_section {
public:
	n_tar(){
                mass = 0.939565;
                PID = 2112;
                mag_mom = -1.913;
        }
};

std::vector<std::string> valid_targets = {"p_tar","n_tar"};

#endif // TARGETS_H

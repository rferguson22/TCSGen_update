#ifndef TARGETS_H
#define TARGETS_H

#include "cross_sections.h"  // Include the base class
#include <vector>
#include <string>
#include <iostream>

class ClassRegistry {
public:
	static void registerClass(const std::string& className) {
		 getClassNames().push_back(className);
	 }
    
	static const std::vector<std::string>& getClassNames() {
		return getClassNamesInternal();
	 }
private:
	static std::vector<std::string>& getClassNamesInternal() {
		static std::vector<std::string> classNames;
		 return classNames;
	}
};

class p_tar : public cross_sections {
public:
	mass = 0.938;
	PID = 2212;
	mag_mom = 1.410;

	p_tar() {
		ClassRegistry::registerClass(className());
	}

	virtual std::string className() const override { return "p_tar"; }
};

class n_tar : public cross_sections {
public:
	mass = 0.939;
	PID = 2112;
	mag_mom = -1.913;

	n_tar() {
		ClassRegistry::registerClass(className())
	}

	virtual std::string className() const override { return "n_tar"; }
};


#endif // TARGETS_H

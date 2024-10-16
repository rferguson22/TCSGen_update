#ifndef MODELS_H
#define MODELS_H

#include "Targets.h"
#include <memory>
#include <vector>
#include <string>


class model_proton : public p_tar {
public:
	model_proton():p_tar(){}

	static const std::string type_name;
	std::string get_type() const override {
		return type_name;
	}
	
	double f1(double t)override{
		return (1/((1-t/0.71)*(1-t/0.71)))*
			(1/(1-t/(4*mass*mass)))*
			(1-mag_mom*t/(4*mass*mass));
	}

	double f2(double t)override{
		return (1/((1-t/0.71)*(1-t/0.71)))*
			(1/(1-t/(4*mass*mass)))*
			(mag_mom-1);
	}
};

const std::string model_proton::type_name = "model_proton";



class model_proton_test : public p_tar {
public:
	model_proton_test():p_tar(){}	
	
	static const std::string type_name;
        std::string get_type() const override {
                return type_name;
        }

        double f1(double t)override{
                return (1/((1-t/0.9)*(1-t/0.9)))*
                        (1/(1-t/(4*mass*mass)))*
                        (1-mag_mom*t/(4*mass*mass));
        }

        double f2(double t)override{
                return (1/((1-t/0.9)*(1-t/0.9)))*
                        (1/(1-t/(4*mass*mass)))*
                        (mag_mom-1);
        }
};

const std::string model_proton_test::type_name = "model_proton_test";





class model_neutron : public n_tar {
public:
	model_neutron():n_tar(){}
	
	static const std::string type_name;
        std::string get_type() const override {
                return type_name;
        }

        double f1(double t)override{
                return (1/((1-t/0.71)*(1-t/0.71)))*
			(1/(1-t/(4*mass*mass)))*
			(-mag_mom*t/(4*mass*mass));
        }

        double f2(double t)override{
                return (1./((1-t/0.71)))*
			(1/(1-t/(4.*mass*mass)))*
			mag_mom;                        
        }
};

const std::string model_neutron::type_name = "model_neutron";





class model_neutron_test : public n_tar {
public:
	model_neutron_test():n_tar(){}
	
	static const std::string type_name;
        std::string get_type() const override {
                return type_name;
        }

        double f1(double t)override{
                return (1/((1-t/0.9)*(1-t/0.9)))*
                        (1/(1-t/(4*mass*mass)))*
                        (-mag_mom*t/(4*mass*mass));
        }

        double f2(double t)override{
                return (1./((1-t/0.9)))*          
                        (1/(1-t/(4.*mass*mass)))*
                        mag_mom;
        }
};

const std::string model_neutron_test::type_name = "model_neutron_test";


std::map<std::string, std::vector<std::shared_ptr<cross_section>>> target_models = {
	{"p_tar", {std::make_shared<model_proton>(),
		   std::make_shared<model_proton_test>()}},
	{"n_tar", {std::make_shared<model_neutron>(),
		   std::make_shared<model_neutron_test>()}},
};

std::shared_ptr<cross_section>create_model(const std::string& model_type) {

	if (model_type == "model_proton"){
		return std::make_shared<model_proton>();
	}
	else if (model_type == "model_proton_test"){
		return std::make_shared<model_proton_test>();
	}
	else if (model_type == "model_neutron"){
                return std::make_shared<model_neutron>();
        }
	else if (model_type == "model_neutron_test"){
                return std::make_shared<model_neutron_test>();
        }
	else{
		return nullptr;
	}
}



#endif // MODELS_H

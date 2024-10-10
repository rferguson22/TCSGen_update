#ifndef MODELS_H
#define MODELS_H

#include "Targets.h"
#include <memory>
#include <vector>
#include <string>


class model_proton : public p_tar {
public:
	model_proton():p_tar(){}
	
	static const string type_name;

	string get_type() const override {
		return type_name;
	}

	double f1(double t) const override{
		return (1/((1-t/0.71)*(1-t/0.71)))*
			(1/(1-t/(4*mass*mass)))*
			(1-mag_mom*t/(4*mass*mass));
	}

	double f2(double t) const override{
		return (1/((1-t/0.71)*(1-t/0.71)))*
			(1/(1-t/(4*mass*mass)))*
			(mag_mom-1);
	}
};

const string model_proton::type_name = "model_proton";



class model_proton_test : public p_tar {
public:
	model_proton_test():p_tar(){}	

        static const string type_name;

        string get_type() const override {
                return type_name;
        }

        double f1(double t) const override{
                return (1/((1-t/0.9)*(1-t/0.9)))*
                        (1/(1-t/(4*mass*mass)))*
                        (1-mag_mom*t/(4*mass*mass));
        }

        double f2(double t) const override{
                return (1/((1-t/0.9)*(1-t/0.9)))*
                        (1/(1-t/(4*mass*mass)))*
                        (mag_mom-1);
        }
};

const string model_proton_test::type_name = "model_proton_test";





class model_neutron : public n_tar {
public:
	model_neutron():n_tar(){}

        static const string type_name;

        string get_type() const override {
                return type_name;
        }

        double f1(double t) const override{
                return (1/((1-t/0.71)*(1-t/0.71)))*
			(1/(1-t/(4*mass*mass)))*
			(-mag_mom*t/(4*mass*mass));
        }

        double f2(double t) const override{
                return (1./((1-t/0.71)))*
			(1/(1-t/(4.*mass*mass)))*
			mag_mom;                        
        }
};

const string model_neutron::type_name = "model_neutron";





class model_neutron_test : public n_tar {
public:
	model_neutron_test():n_tar(){}

        static const string type_name;

        string get_type() const override {
                return type_name;
        }

        double f1(double t) const override{
                return (1/((1-t/0.9)*(1-t/0.9)))*
                        (1/(1-t/(4*mass*mass)))*
                        (-mag_mom*t/(4*mass*mass));
        }

        double f2(double t) const override{
                return (1./((1-t/0.9)))*          
                        (1/(1-t/(4.*mass*mass)))*
                        mag_mom;
        }
};

const string model_neutron_test::type_name = "model_neutron_test";




shared_ptr<cross_section> cross_target::create_target(const string& c_sec_type,const string& model_type) {
	shared_ptr<cross_section> target;

	if (model_type == "model_proton"){
		return make_shared<model_proton>();
	}
	else if (model_type == "model_proton_test"){
		return make_shared<model_proton_test>();
	}
	else if (model_type == "model_neutron"){
                return make_shared<model_neutron>();
        }
	else if (model_type == "model_neutron_test"){
                return make_shared<model_neutron_test>();
        }
	else{
		return nullptr;
	}

	return target;
}


vector<string> get_valid_models(){
	return {model_proton::type_name, 
		model_proton_test::type_name,
		model_neutron::type_name,
		model_neutron_test::type_name};
}


#endif // MODELS_H

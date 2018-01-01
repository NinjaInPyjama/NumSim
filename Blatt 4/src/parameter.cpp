#include "parameter.hpp"
#include <iostream>
#include <random>
#include <chrono>

using namespace std;
/// Constructs a new Parameter set with default values
// Driven Cavity parameters; see exercise sheet 1
Parameter::Parameter() {
        
    // load parameters from file
    Load("default.param");
    // Load("actual.param");
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution(1500.0,1000.0/6.0);
    
    _re = distribution(generator);
    
    std::cout << "Re= " << _re << std::endl;

    
}

/// Loads the parameter values from a file
void Parameter::Load(const char * file) {

	FILE* handle = fopen(file, "r");
	double inval;
	char name[20];
	while (!feof(handle)) {
		if (!fscanf(handle, "%s = %lf\n", name, &inval)) continue;
		if (strcmp(name, "re") == 0) _re = inval;
		else if (strcmp(name, "omg") == 0) _omega = inval;
		else if (strcmp(name, "alpha") == 0) _alpha = inval;
		else if (strcmp(name, "dt") == 0) _dt = inval;
		else if (strcmp(name, "tend") == 0) _tend = inval;
		else if (strcmp(name, "iter") == 0) _itermax = inval;
		else if (strcmp(name, "eps") == 0) _eps = inval;
		else if (strcmp(name, "tau") == 0) _tau = inval;
		else printf("Unknown parameter %s\n", name);
	}
	fclose(handle);
}

const real_t &Parameter::Re() const{
    return _re;
}
const real_t &Parameter::Omega() const{
    return _omega;
}
const real_t &Parameter::Alpha() const{
    return _alpha;
}
const real_t &Parameter::Dt() const{
    return _dt;
}
const real_t &Parameter::Tend() const{
    return _tend;
}
const index_t &Parameter::IterMax() const{
    return _itermax;
}
const real_t &Parameter::Eps() const{
    return _eps;
}
const real_t &Parameter::Tau() const{
    return _tau;
}

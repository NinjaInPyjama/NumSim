#include "parameter.hpp"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;
/// Constructs a new Parameter set with default values
// Driven Cavity parameters; see exercise sheet 1
Parameter::Parameter() {
    
    //parameter aus example 1.3 
    // dt, tend geraten
    
    _re = 1000;
    _omega = 1.7;
    _alpha = 0.9;
    _dt = 0.05;
    _tend = 50;
    _eps = 0.001;
    _tau = 0.5;
    _itermax = 100;
    
    
    // eigentlich mit if file exits
    this->Load("default.param");
    
    //this->Load("actual.param");
   
    
}

/// Loads the parameter values from a file
void Parameter::Load(const char * file) {
    
    ifstream fin(file);
    real_t a;
    string name;
    string gleich;
    while (fin >> name >> gleich >> a){
        //cout << name << gleich << a;
        if (name.compare("eps")){
			_eps = a;
        }
        else if (name.compare("alpha")){
			_alpha = a;
        }
		else if (name.compare("omg")){
			_omega = a;
        }
		else if (name.compare("re")){
			_re = a;
        }
        else if (name.compare("dt")){
			_dt = a;	
        }
        else if (name.compare("tend")){
			_tend = a;	
        }
        else if (name.compare("iter")){
			_itermax = index_t(a);
        }
        else if (name.compare("tau")){
			_tau = a;
        }
		else {
			cout << "unknown parameter name" << endl;
        }
    }
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

#include "compute.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include <math.h>
#include "solver.hpp"
#include "parameter.hpp"

/// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry * geom, const Parameter * param) {
    
    _geom = geom;
    _param = param;
    
    Grid * _tmp = new Grid(_geom);
    
    Grid * _F = new Grid(_geom, multi_real_t(1.0, 0.5));
    
    Grid * _G = new Grid(_geom, multi_real_t(0.5, 1.0));
    
    Grid * _p = new Grid(_geom, multi_real_t(0.5, 0.5));
    
    Grid * _v = new Grid(_geom, multi_real_t(0.5, 1.0));
    
    Grid * _u = new Grid(_geom, multi_real_t(1.0, 0.5));
    
    Solver * _solver = new SOR(_geom, _param->Omega());
    
    
    real_t _t = 0.0;

    real_t _dtlimit = _param->Dt();

    real_t _epslimit = _param->Eps();



  
    
}

/// Deletes all grids
Compute::~Compute() {

}


/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {

}

/// Returns the simulated time in total
const real_t & Compute::GetTime() const {
	return _t;
}


/// Returns the pointer to U
const Grid * Compute::GetU() const {
	return _u;
}

/// Returns the pointer to V
const Grid * Compute::GetV() const {
	return _v;
}

/// Returns the pointer to P
const Grid * Compute::GetP() const {
	return _p;
}

/// Returns the pointer to RHS
const Grid * Compute::GetRHS() const {
	return _rhs;
}


/// Computes and returns the absolute velocity
const Grid * Compute::GetVelocity() {
    InteriorIterator iit = InteriorIterator(_geom);
    Grid * abs_vel = new Grid(_geom);
    multi_real_t mitte_der_zelle = multi_real_t(0.5, 0.5);
    real_t v_interpol = 0.0;
    real_t u_interpol = 0.0;
    for(iit.First();iit.Valid();iit.Next()){
        
        v_interpol = _v->Interpolate(multi_real_t(real_t(iit.Pos()[0]) + mitte_der_zelle[0] , real_t(iit.Pos()[1]) + mitte_der_zelle[1] ));
        u_interpol = _u->Interpolate(multi_real_t(real_t(iit.Pos()[0]) + mitte_der_zelle[0] , real_t(iit.Pos()[1]) + mitte_der_zelle[1] ));
        abs_vel->Cell(iit) = sqrt(v_interpol*v_interpol + u_interpol*u_interpol);
    }
    return abs_vel;
}

/// Computes and returns the vorticity
const Grid * Compute::GetVorticity() {
        
    InteriorIterator iit = InteriorIterator(_geom);
    
    Grid * du_dy = new Grid(_geom);
    Grid * dv_dx = new Grid(_geom);
    for(iit.First();iit.Valid();iit.Next()){
        du_dy->Cell(iit) = _u->dy_central(iit);
        dv_dx->Cell(iit) = _v->dx_central(iit);
    }
    
    Grid * vorticity = new Grid(_geom);
    multi_real_t mitte_der_zelle = multi_real_t(0.5, 0.5);
    for(iit.First();iit.Valid();iit.Next()){
        
        vorticity->Cell(iit) = dv_dx->Interpolate(multi_real_t(real_t(iit.Pos()[0]) + mitte_der_zelle[0] , real_t(iit.Pos()[1]) + mitte_der_zelle[1] )) - du_dy->Interpolate(multi_real_t(real_t(iit.Pos()[0]) + mitte_der_zelle[0] , real_t(iit.Pos()[1]) + mitte_der_zelle[1] ));
        
    }
    return vorticity;
    
}

/// Computes and returns the stream line values
const Grid * Compute::GetStream() {
	return nullptr;
}


/// Compute the new velocites u,v
void Compute::NewVelocities(const real_t & dt) {
    
    InteriorIterator iit = InteriorIterator(_geom);
    for(iit.First();iit.Valid();iit.Next()){
        
        _u->Cell(iit) = _u->Cell(iit) - dt*(_p->dx_r(iit));
        _v->Cell(iit) = _v->Cell(iit) - dt*(_p->dy_l(iit));
        
    }
    
}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);
    for(iit.First();iit.Valid();iit.Next()){
        
        _F->Cell(iit) = _u->Cell(iit) + dt*(( _u->dxx(iit) + _u->dyy(iit) )/_param->Re() - _u->DC_udu_x(iit, _param->Alpha()) - _u->DC_vdu_y(iit, _param->Alpha(), _v));
        _G->Cell(iit) = _v->Cell(iit) + dt*(( _v->dxx(iit) + _v->dyy(iit) )/_param->Re() - _v->DC_vdv_y(iit, _param->Alpha()) - _v->DC_udv_x(iit, _param->Alpha(), _u));
        
    }
}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t & dt) {
    
    InteriorIterator iit = InteriorIterator(_geom);
    for(iit.First();iit.Valid();iit.Next()){
        
        //_rhs->Cell(iit) = (_p->dxx(iit) + _p->dyy(iit)); 
        _rhs->Cell(iit) = (_F->dx_l(iit) + _G->dy_r(iit))/dt ;
        
    }
    
}

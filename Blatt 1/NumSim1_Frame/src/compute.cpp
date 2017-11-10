#include "compute.hpp"

using namespace std;

/// Creates a compute instance with given geometry and parameter
Compute::Compute(const Geometry * geom, const Parameter * param) {
    
    _geom = geom;
    _param = param;
    
    _tmp = new Grid(_geom);
    
    _F = new Grid(_geom, multi_real_t(1.0, 0.5));
	_F->Initialize(0.0);
    _G = new Grid(_geom, multi_real_t(0.5, 1.0));
	_G->Initialize(0.0);

    _p = new Grid(_geom, multi_real_t(0.5, 0.5));
    _p->Initialize(_geom->Pressure());
    _v = new Grid(_geom, multi_real_t(0.5, 1.0));
    _v->Initialize(_geom->Velocity()[1]);
    _u = new Grid(_geom, multi_real_t(1.0, 0.5));
	_u->Initialize(_geom->Velocity()[0]);
	_geom->Update_U(_u);

	_rhs = new Grid(_geom, multi_real_t(0.5, 0.5));
	_rhs->Initialize(0.0);

	_solver = new SOR(_geom); //, _param->Omega());
        
    _t = 0.0;
    _dtlimit = _param->Dt();
    _epslimit = _param->Eps();
}

/// Deletes all grids
Compute::~Compute() {

}


/// Execute one time step of the fluid simulation (with or without debug info)
// @ param printInfo print information about current solver state (residual
// etc.)
void Compute::TimeStep(bool printInfo) {
    // TODO: see script p. 23
    
	// boundary_val(...)
	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);

    // compute_fg(...)
    MomentumEqu(_dtlimit);
	_geom->Update_U(_F);
	_geom->Update_V(_G);
	//_F->print();

    // compute_rhs(...)
    RHS(_dtlimit);

    // solver iteration
    index_t itermax = _param->IterMax();
    index_t it = 0;
    real_t res = 0.0;
    do {
        it++;
        res = _solver->Cycle(_p, _rhs);
		if(printInfo) std::cout << "Residual at iteration " << it << ": " << res << std::endl;
    } while(it<itermax && res>_epslimit);
    if(printInfo) std::cout << "Solver stopped at iteration " << it << " with residual: " << res << std::endl;

	//_p->print();

    // compute_uv(...)
    NewVelocities(_dtlimit);

	// boundary_val(...)
	_geom->Update_U(_u);
	_geom->Update_V(_v);
	_geom->Update_P(_p);
	
    _t += _dtlimit;
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
	multi_real_t cell_center = multi_real_t(0.5*_geom->Mesh()[0], 0.5*_geom->Mesh()[1]);
	real_t u_ip = 0.0; // storage for interpolated u to center of cells
	real_t v_ip = 0.0; // storage for interpolated v to center of cells


    for(iit.First(); iit.Valid(); iit.Next()){
        // Interpolating the velocities to center of cells
        multi_index_t cell_pos = iit.Pos();
		v_ip = (_v->Cell(iit.Down()) + _v->Cell(iit)) / 2.0; // _v->Interpolate(multi_real_t((real_t)cell_pos[0] + cell_center[0], (real_t)cell_pos[1] + cell_center[1]));
		u_ip = (_u->Cell(iit.Left()) + _u->Cell(iit)) / 2.0; // _u->Interpolate(multi_real_t((real_t)cell_pos[0] + cell_center[0], (real_t)cell_pos[1] + cell_center[1]));
        abs_vel->Cell(iit) = sqrt(v_ip*v_ip + u_ip*u_ip);
    }

	BoundaryIterator bit = BoundaryIterator(_geom);
	bit.SetBoundary(0);
	abs_vel->Cell(bit.Left()) = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		abs_vel->Cell(bit) = 0.0;
	}
	// Iteration over right boundary
	bit.SetBoundary(1);
	abs_vel->Cell(bit.Top()) = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		abs_vel->Cell(bit) = 0.0;
	}
	// Iteration over lower boundary
	bit.SetBoundary(2);
	abs_vel->Cell(bit.Right()) = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		abs_vel->Cell(bit) = 0.0;
	}

	// Iteration over left boundary
	bit.SetBoundary(3);
	abs_vel->Cell(bit.Down()) = 0;
	for (bit.First(); bit.Valid(); bit.Next()) {
		abs_vel->Cell(bit) = 0.0;
	}
    return abs_vel;
}

/// Computes and returns the vorticity
const Grid * Compute::GetVorticity() {
    InteriorIterator iit = InteriorIterator(_geom);
    Grid * vort = new Grid(_geom);
    multi_real_t cell_center = multi_real_t(0.5, 0.5);

    // Creating grids of derivatives of u (in y-dim) and v (in x-dim) (for interpolation reasons)
    Grid * du_dy = new Grid(_geom);
    Grid * dv_dx = new Grid(_geom);
    for(iit.First();iit.Valid();iit.Next()){
        du_dy->Cell(iit) = _u->dy_central(iit);
        dv_dx->Cell(iit) = _v->dx_central(iit);
    }
    
    for(iit.First();iit.Valid();iit.Next()){ 
        // Calculating vorticity by dv/dx - du/dy
        multi_index_t cell_pos = iit.Pos();
		vort->Cell(iit) = 0.0; //  dv_dx->Interpolate(multi_real_t((real_t)cell_pos[0] + cell_center[0], (real_t)cell_pos[1] + cell_center[1]))
                                - du_dy->Interpolate(multi_real_t((real_t)cell_pos[0] + cell_center[0], (real_t)cell_pos[1]  + cell_center[1]));
    }
    return vort;
}

/// Computes and returns the stream line values
const Grid * Compute::GetStream() {
    // whatever
	return nullptr;
}


/// Compute the new velocites u,v
void Compute::NewVelocities(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);
    
    for(iit.First();iit.Valid();iit.Next()){
		// see script, p. 18, formular (3.1)
        _u->Cell(iit) = _F->Cell(iit) - dt*(_p->dx_r(iit));
        _v->Cell(iit) = _G->Cell(iit) - dt*(_p->dy_r(iit));
    }
}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);

    for(iit.First();iit.Valid();iit.Next()){
		// see script, p. 18, formular (3.2)
		real_t A = (_u->dxx(iit) + _u->dyy(iit)) / _param->Re() - _u->DC_du2_x(iit, _param->Alpha()) - _u->DC_duv_y(iit, _param->Alpha(), _v);
		real_t B = (_v->dxx(iit) + _v->dyy(iit)) / _param->Re() - _v->DC_dv2_y(iit, _param->Alpha()) - _v->DC_duv_x(iit, _param->Alpha(), _u);
		_F->Cell(iit) = _u->Cell(iit) + dt*A;
		_G->Cell(iit) = _v->Cell(iit) + dt*B;
    }
}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t & dt) {
    InteriorIterator iit = InteriorIterator(_geom);

    for(iit.First();iit.Valid();iit.Next()){
		// see script, p. 19, formular (3.3)
        _rhs->Cell(iit) = (_F->dx_l(iit) + _G->dy_l(iit))/dt;
    }
}